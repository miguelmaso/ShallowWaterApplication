# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MultiGridRefiningProcess(model, settings["Parameters"])

class MultiGridRefiningProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)
        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"            : "MainModelPart",
            "maximum_number_of_subscales"     : 4,
            "current_subscale_index"          : 0,
            "refining_condition_parameters"   : {
                "kratos_module"                   : "",
                "python_module"                   : "",
                "Parameters"                      : {}
            },
            "variables_to_apply_fixity"       : [],
            "variables_to_set_at_interface"   : [],
            "variables_to_update_at_coarse"   : [],
            "advanced_configuration"          : {
                "echo_level"                      : 0,
                "number_of_divisions_at_subscale" : 2,
                "subscale_interface_base_name"    : "refined_interface",
                "subscale_boundary_condition"     : "LineCondition2D2N"
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)

        self.model = model

        self.number_of_divisions_at_subscale = self.settings['advanced_configuration']['number_of_divisions_at_subscale'].GetInt()
        self.number_of_substeps = 2**self.settings['advanced_configuration']['number_of_divisions_at_subscale'].GetInt()

        current_subscale_index = self.settings['current_subscale_index'].GetInt()

        # Get the coarse and refined model part names
        model_part_name = self.settings['model_part_name'].GetString()
        self.coarse_model_part_name = model_part_name + '_' + str(current_subscale_index)
        self.refined_model_part_name = model_part_name + '_' + str(current_subscale_index + 1)

        self.coarse_model_part = self.model.GetModelPart(self.coarse_model_part_name)
        self.refined_model_part = self.model.CreateModelPart(self.refined_model_part_name)
        self.refined_model_part.ProcessInfo[MeshingApplication.SUBSCALE_INDEX] = current_subscale_index + 1

        self.variables_to_apply_fixity = GenerateVariableListFromInput(self.settings["variables_to_apply_fixity"])
        self.variables_to_set_at_interface = GenerateVariableListFromInput(self.settings["variables_to_set_at_interface"])
        self.variables_to_update_at_coarse = GenerateVariableListFromInput(self.settings["variables_to_update_at_coarse"])

        kratos_module_name = self.settings["refining_condition_parameters"]["kratos_module"].GetString()
        python_module_name = self.settings["refining_condition_parameters"]["python_module"].GetString()
        if kratos_module_name:
            full_module_name = kratos_module_name + "." + python_module_name
        else:
            full_module_name = python_module_name
        python_module = __import__(full_module_name, fromlist=[python_module_name])
        self.settings["refining_condition_parameters"]["Parameters"]["model_part_name"].SetString(self.coarse_model_part_name)
        self.refining_condition_process = python_module.Factory(self.settings["refining_condition_parameters"], self.model)

    def ExecuteInitialize(self):
        self._InitializeRefinedModelPart()
        self.subscales_utility = MeshingApplication.MultiscaleRefiningProcess(
            self.coarse_model_part,
            self.refined_model_part,
            self.model.CreateModelPart("to_delete"), # this was the visualization model part. TODO: remove
            self.settings["advanced_configuration"])

    def ExecuteBeforeSolutionLoop(self):
        self._EvaluateCondition()
        self._ExecuteRefinement()
        self._ExecuteCoarsening()

    def ExecuteInitializeSolutionStep(self):
        self._TransferSubstepToRefinedInterface()
        self._ApplyFixityAtInterface(True)

    def ExecuteFinalizeSolutionStep(self):
        self._ApplyFixityAtInterface(False)
        self._UpdateVariablesAtCoarseModelPart()

    def _InitializeRefinedModelPart(self):
        buffer_size = self.coarse_model_part.GetBufferSize()
        self.refined_model_part.SetBufferSize(buffer_size)
        MeshingApplication.MultiscaleRefiningProcess.InitializeRefinedModelPart(self.coarse_model_part, self.refined_model_part)

    def _ExecuteRefinement(self):
        self.subscales_utility.ExecuteRefinement()

    def _ExecuteCoarsening(self):
        self.subscales_utility.ExecuteCoarsening()

    def _EvaluateCondition(self):
        self.refining_condition_process.Execute()

    def _ApplyFixityAtInterface(self, state):
        for variable in self.variables_to_apply_fixity:
            self.subscales_utility.FixRefinedInterface(variable, state)

    def _TransferSubstepToRefinedInterface(self):
        substep_fraction = self.refined_model_part.ProcessInfo[KratosMultiphysics.STEP] / self.number_of_substeps
        for variable in self.variables_to_set_at_interface:
            self.subscales_utility.TransferSubstepToRefinedInterface(variable, substep_fraction)

    def _UpdateVariablesAtCoarseModelPart(self):
        for variable in self.variables_to_update_at_coarse:
            self.subscales_utility.TransferLastStepToCoarseModelPart(variable)
