# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.MeshingApplication as MSH

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.wave_solver import WaveSolver
from KratosMultiphysics.MeshingApplication.multiscale_refining_process import MultiscaleRefiningProcess

def CreateSolver(model, custom_settings):
    return MultigridSolver(model, custom_settings)

class MultigridSolver(WaveSolver):

    def __init__(self, model, settings):
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model = model      # TODO: inherit from PythonSolver and use super
        self.settings = settings
        self.echo_level = self.settings["echo_level"].GetInt()

        ## Set the element and condition names for the replace settings
        self.element_name, self.condition_name, self.min_buffer_size = self._GetFormulationSettings()

        # Initialize the multigrid process. It creates the model part
        self.multigrid = MultiscaleRefiningProcess(model, settings["multigrid_settings"])
        self.main_model_part = self.multigrid.GetRefinedModelPart()

        self._SetProcessInfo()

    def ImportModelPart(self):
        if self.main_model_part.ProcessInfo[MSH.SUBSCALE_INDEX] == 0:
            # Default implementation in the base class
            self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if self.main_model_part.ProcessInfo[MSH.SUBSCALE_INDEX] == 0:
            super().PrepareModelPart()
        self.multigrid.PrepareModelPart() # It creates the cpp utility instance

    def AdvanceInTime(self, current_time):
        divisions = 2**(self.GetComputingModelPart().GetValue(MSH.SUBSCALE_INDEX) * self.multigrid.number_of_divisions_at_subscale)
        current_time += self._GetEstimateDeltaTimeUtility().Execute() / divisions

        self.GetComputingModelPart().CloneTimeStep(current_time)
        self.GetComputingModelPart().ProcessInfo[KM.STEP] += 1

        self.multigrid.ExecuteInitializeSolutionStep()

        return current_time

    def GetComputingModelPart(self):
        return self.multigrid.GetRefinedModelPart()

    def GetStepLabel(self):
        if self.GetComputingModelPart().GetValue(MSH.SUBSCALE_INDEX) == 0:
            return self.GetComputingModelPart().ProcessInfo[KM.STEP]
        else:
            level = self.GetComputingModelPart().GetValue(MSH.SUBSCALE_INDEX)
            coarse_step = self.multigrid.GetCoarseModelPart().ProcessInfo[KM.STEP]
            refined_step = self.GetComputingModelPart().ProcessInfo[KM.STEP]
            return f'[{level}] {coarse_step}.{refined_step}'

    def InitializeSolutionStep(self):
        if self.GetComputingModelPart().NumberOfElements() != 0:
            super().InitializeSolutionStep()

    def Predict(self):
        if self.GetComputingModelPart().NumberOfElements() != 0:
            super().Predict()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            if self.GetComputingModelPart().NumberOfElements() != 0:
                is_converged = self._GetSolutionStrategy().SolveSolutionStep()
                return is_converged

    def FinalizeSolutionStep(self):
        if self.GetComputingModelPart().NumberOfElements() != 0:
            super().FinalizeSolutionStep()

    def Finalize(self):
        self.multigrid.ExecuteFinalize()

    def GetDefaultParameters(self):
        default_parameters = KM.Parameters("""
        {
            "multigrid_settings" : {}
        }
        """)
        default_parameters.AddMissingParameters(super().GetDefaultParameters())
        return default_parameters
