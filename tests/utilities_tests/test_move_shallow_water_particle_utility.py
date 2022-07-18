import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.process_factory import KratosProcessFactory as ProcessFactory

from math import sin, pi

class TestMoveShallowWaterParticleUtility(KratosUnittest.TestCase):
    """Wrapper for testing the move shallow water particles utility"""

    def testPureConvection(self):
        model = KM.Model()
        settings = KM.Parameters('''{
            "particles_settings" : {
                "convection_scalar_variable"  : "HEIGHT",
                "convection_vector_variable"  : "VELOCITY"
            },
            "processes_list" : [{
                "kratos_module"   : "KratosMultiphysics",
                "python_module"   : "from_json_check_result_process",
                "Parameters"      : {
                    "model_part_name"  : "model_part",
                    "check_variables"  : ["HEIGHT"],
                    "input_file_name"  : "utilities_tests/test_pure_convection_reference.json",
                    "time_frequency"   : 0.4
                }
            }],
            "output_processes_list" : []
        }''')
        test = AuxiliaryAnalysis(model, settings)

        for node in test.model_part.Nodes:
            node.SetSolutionStepValue(KM.VELOCITY, self._Vortex(node.X, node.Y))
            if ((node.X - 0.4)**2 + (node.Y - 0.4)**2) < 0.01:
                node.SetSolutionStepValue(SW.HEIGHT, 0, 1.0)
                node.SetSolutionStepValue(SW.HEIGHT, 1, 1.0)
            else:
                node.SetSolutionStepValue(SW.HEIGHT, 0, 0.8)
                node.SetSolutionStepValue(SW.HEIGHT, 1, 0.8)

        test.Initialize()
        test.RunSolutionLoop()
    
    def testReseed(self):
        model = KM.Model()
        settings = KM.Parameters('''{
            "particles_settings" : {
                "convection_scalar_variable"  : "HEIGHT",
                "convection_vector_variable"  : "VELOCITY"
            },
            "processes_list" : [{
                "kratos_module"   : "KratosMultiphysics",
                "python_module"   : "from_json_check_result_process",
                "Parameters"      : {
                    "model_part_name"  : "model_part",
                    "check_variables"  : ["HEIGHT"],
                    "input_file_name"  : "utilities_tests/test_reseed_reference.json",
                    "time_frequency"   : 0.4
                }
            }],
            "output_processes_list" : []
        }''')
        test = AuxiliaryAnalysis(model, settings)

        for node in test.model_part.Nodes:
            node.SetSolutionStepValue(KM.VELOCITY, [1.0, 0.0, 0.0])
            if node.X < 0.01:
                node.SetSolutionStepValue(SW.HEIGHT, 1.0)
            else:
                node.SetSolutionStepValue(SW.HEIGHT, 0.0)

        test.Initialize()
        test.RunSolutionLoop()

    @staticmethod
    def _Vortex(x, y):
        L = 1.0
        x0 = 0.5 * L
        y0 = 0.5 * L
        return [0.2 * sin(pi*x/L)*sin(2*pi*y/L), -0.2 * sin(2*pi*x/L)*sin(pi*y/L), 0.0]

    auxiliary_settings = KM.Parameters('''{
        "processes_list_not_to_execute" : [{
            "kratos_module"   : "KratosMultiphysics",
            "python_module"   : "json_output_process",
            "Parameters"      : {
                "model_part_name"  : "model_part",
                "output_variables" : ["HEIGHT"],
                "output_file_name" : "utilities_tests/test_XXX_reference.json",
                "time_frequency"   : 0.4
            }
        },{
            "kratos_module"   : "KratosMultiphysics",
            "python_module"   : "from_json_check_result_process",
            "Parameters"      : {
                "model_part_name"  : "model_part",
                "check_variables"  : ["HEIGHT"],
                "input_file_name"  : "utilities_tests/test_XXX_reference.json",
                "time_frequency"   : 0.4
            }
        },{
            "kratos_module"        : "KratosMultiphysics",
            "python_module"        : "gid_output_process",
            "Parameters"           : {
                "model_part_name"        : "model_part",
                "output_name"            : "utilities_tests/test_particles",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "output_control_type"   : "step",
                        "output_interval"       : 1,
                        "nodal_results"         : ["HEIGHT","VELOCITY"]
                    }
                }
            }
        }]
    }''')

class AuxiliaryAnalysis:
    """Auxiliary class to execute the particles utility."""

    def __init__(self, model, parameters):
        self.model = model
        self.parameters = parameters
        self.model_part = self.model.CreateModelPart("model_part")
        self.model_part.SetBufferSize(2)

        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)
        self.model_part.AddNodalSolutionStepVariable(SW.HEIGHT)
        self.model_part.AddNodalSolutionStepVariable(SW.DELTA_SCALAR)
        self.model_part.AddNodalSolutionStepVariable(SW.DELTA_VECTOR)
        self.model_part.AddNodalSolutionStepVariable(SW.PROJECTED_SCALAR)
        self.model_part.AddNodalSolutionStepVariable(SW.PROJECTED_VECTOR)
        self.model_part.AddNodalSolutionStepVariable(KM.INTEGRATION_WEIGHT)

        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 2
        self.model_part.ProcessInfo[KM.GRAVITY_Z] = 9.81

        KM.ModelPartIO("processes_tests/model_part").ReadModelPart(self.model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.model_part)

    def Initialize(self):

        self.step = 0
        self.time = 0.0
        self.end_time = 0.5
        self.time_step = 0.1

        self.model_part.ProcessInfo[KM.STEP] = self.step
        self.model_part.ProcessInfo[KM.TIME] = self.time

        factory = ProcessFactory(self.model)
        self.processes = factory.ConstructListOfProcesses(self.parameters["processes_list"])
        self.output_processes = factory.ConstructListOfProcesses(self.parameters["output_processes_list"])
        self.processes += self.output_processes

        KM.FindGlobalNodalNeighboursProcess(self.model_part).Execute()
        KM.GenericFindElementalNeighboursProcess(self.model_part).ExecuteInitialize()

        self.particles = SW.MoveShallowWaterParticleUtility(self.model_part, self.parameters["particles_settings"])
        self.particles.MountBin()
        self.pre_minimum_num_of_particles = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        self.post_minimum_num_of_particles = 2 * self.model_part.ProcessInfo[KM.DOMAIN_SIZE]

        for process in self.processes:
            process.ExecuteInitialize()
            process.ExecuteBeforeSolutionLoop()

    def _Solve(self):
        pass

    def _MoveParticles(self):
        self.particles.CalculateVelOverElemSize()
        self.particles.MoveParticles()
        self.particles.PreReseed(self.pre_minimum_num_of_particles)
        self.particles.TransferLagrangianToEulerian()
        KM.VariableUtils().CopyVariable(SW.PROJECTED_SCALAR, SW.HEIGHT, self.model_part.Nodes)
        self.particles.CopyScalarVarToPreviousTimeStep(SW.HEIGHT, self.model_part.Nodes)
        self.particles.ResetBoundaryConditions()

    def _UpdateParticles(self):
        self.particles.CalculateDeltaVariables()
        self.particles.CorrectParticlesWithoutMovingUsingDeltaVariables()
        self.particles.PostReseed(self.post_minimum_num_of_particles)

    def RunSolutionLoop(self):
        while self.time < self.end_time:
            self.step += 1
            self.time += self.time_step

            self.model_part.CloneTimeStep(self.time)
            self.model_part.ProcessInfo[KM.STEP] = self.step

            for process in self.processes:
                process.ExecuteInitializeSolutionStep()

            self._MoveParticles()
            self._Solve()
            self._UpdateParticles()

            for process in self.processes:
                process.ExecuteFinalizeSolutionStep()

            for process in self.output_processes:
                process.PrintOutput()

if __name__ == '__main__':
    KratosUnittest.main()
