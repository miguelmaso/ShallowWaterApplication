# Importing Kratos
import KratosMultiphysics as KM

from KratosMultiphysics.ShallowWaterApplication.shallow_water_analysis import ShallowWaterAnalysis
from KratosMultiphysics.ShallowWaterApplication.legacy.multigrid_refining_process import MultiGridRefiningProcess

class ShallowWaterMultigridAnalysis(ShallowWaterAnalysis):
    ''' Main script for shallow water simulations '''

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

        self.current_subscale = 0
        self.maximum_subgrids = 0
        if self.project_parameters.Has("multigrid_settings"):
            if not self.project_parameters["multigrid_settings"].Has("maximum_number_of_subscales"):
                raise Exception("Multigrid settings defined but the 'maximum_number_of_subscales' is not present")
            if not self.project_parameters["multigrid_settings"].Has("current_subscale_index"):
                raise Exception("Multigrid settings defined but the 'current_subscale_index' is not present")
            self.maximum_subgrids = self.project_parameters["multigrid_settings"]["maximum_number_of_subscales"].GetInt()
            self.current_subscale = self.project_parameters["multigrid_settings"]["current_subscale_index"].GetInt()

            if self.current_subscale == 0:
                name = self._GetSolver().GetComputingModelPart().Name
                model.RenameModelPart(name, name + '_0')

        if self.current_subscale < self.maximum_subgrids:
            # Initialize the multigrid process. It creates the model part
            self.multigrid = MultiGridRefiningProcess(model, self.project_parameters["multigrid_settings"])
            # Create a nested analysis
            sub_parameters = self._CreateSubAnalysisParameters()
            self.sub_analysis = ShallowWaterMultigridAnalysis(self.model, sub_parameters)

        self._UpdateModelPartNamesInParameters()


    def Initialize(self):
        super().Initialize()
        if self.current_subscale < self.maximum_subgrids:
            self._GetMultiGrid().ExecuteInitialize()
            self.sub_analysis.Initialize()


    def RunSolutionLoop(self):
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.ExecuteCurrentStep()


    def ExecuteCurrentStep(self):
        self.InitializeSolutionStep()
        self._GetSolver().Predict()
        self._GetSolver().SolveSolutionStep()
        self.RunNestedLoop()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()


    def RunNestedLoop(self):
        if self.current_subscale < self.maximum_subgrids:
            self._GetMultiGrid().ExecuteBeforeSolutionLoop()
            number_of_substeps = self._GetMultiGrid().number_of_substeps
            time_step = self._GetSolver().GetComputingModelPart().ProcessInfo[KM.DELTA_TIME]
            nested_loop_time = self.time - time_step
            for _ in range(number_of_substeps):
                nested_loop_time += time_step / number_of_substeps
                self.sub_analysis._GetSolver().AdvanceInTime(nested_loop_time)
                self._GetMultiGrid().ExecuteInitializeSolutionStep()
                self.sub_analysis.ExecuteCurrentStep()
            self._GetMultiGrid().ExecuteFinalizeSolutionStep()


    def _GetMultiGrid(self):
        if hasattr(self, 'multigrid'):
            return self.multigrid
        else:
            return KM.Process()


    def _CreateSubAnalysisParameters(self):
        sub_parameters = self.project_parameters.Clone()
        sub_parameters["multigrid_settings"]["current_subscale_index"].SetInt(self.current_subscale + 1)
        return sub_parameters


    def _UpdateModelPartNamesInParameters(self):
        old_name = self.project_parameters["solver_settings"]["model_part_name"].GetString()
        new_name = self._GetSolver().GetComputingModelPart().Name
        self._UpdateModelPartNameInProcesses(old_name, new_name, self.project_parameters["processes"])
        self._UpdateModelPartNameInProcesses(old_name, new_name, self.project_parameters["output_processes"])
        suffix = new_name[-2:-1]
        self._UpdateOutputNameInProcesses(suffix, self.project_parameters["output_processes"])


    @staticmethod
    def _UpdateModelPartNameInProcesses(old_name, new_name, processes):
        for process_list in processes.values():
            for process in process_list:
                if process.Has("Parameters"):
                    if process["Parameters"].Has("model_part_name"):
                        name = process["Parameters"]["model_part_name"].GetString()
                        name = name.replace(old_name, new_name)
                        process["Parameters"]["model_part_name"].SetString(name)


    @staticmethod
    def _UpdateOutputNameInProcesses(suffix, processes):
        for process_list in processes.values():
            for process in process_list:
                if process.Has("Parameters"):
                    if process["Parameters"].Has("output_name"):
                        name = process["Parameters"]["output_name"].GetString()
                        name = name + suffix
                        process["Parameters"]["output_name"].SetString(name)


    def _GetSimulationName(self):
        return "Shallow Water Multigrid Analysis"
