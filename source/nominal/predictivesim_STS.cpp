/* -------------------------------------------------------------------------- *
 *                     OpenSim:  predictivesim_STS.cpp                        *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                     a                                       *
 * Copyright (c) 2005-2014 Stanford University and the Authors                *
 * Author(s): Carmichael Ong, Ajay Seth                                       *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

//==============================================================================
// This file generates a main() to optimize the set of controls for the torque
// actuators of a model to maximize the distance of a standing long jump with
// soft constraints on landing pose, ligament usage, slipping, and duration. This
// will build an executable that can take the following arguments:
//
// -m <modelfile.osim>: model file to be loaded
// -cf <controlsfile.sto>: initial guess of controls for optimizer
// -opt: flag to run optimization
// -resume: flag to run optimization but resume using resumecmaes.dat file
// -lambda <int>: CMA parameter, set number of samples per outer iteration
// -sigma <double>: CMA parameter, set initial constant of covariance matrix
// -maxIters <int>: optimizer parameter, set max outer iterations
//
// Usage: StandingLongJump_nodal.exe -m AshbyModel_twoConstraints.osim 
//          -cf InitialControllerParameters.sto -opt -lambda 50 -sigma 0.1 -maxIters 2000
//
// The executable assumes that the model file has constraints at the toe and heel
// and that a controller is not present in the model.
//
// The -opt flag is used when an optimization should be run. If the -opt flag is
// not present, then a single forward simulation with analyses will be performed.
//
// Parallelization of the CMA optimizer is done by using MPI. If ENABLE_MPI was
// checked in CMake, then you can call on the executable to use more threads using:
//
// mpiexec -n <numthreads> StandingLongJump_nodal.exe -opt -lambda 50 ...
//==============================================================================

#include <ctime>  // clock(), clock_t, CLOCKS_PER_SEC
#include "EventHandlersAndReporters.h"
#include "ConstraintSlipPenaltyProbe.h"
#include "CoordinateLimitForceProbe.h"
#include "MuscleLikeCoordinateActuator.h"


#include "CMAOptimizer.h"

#ifdef USE_MPI
#include <mpi.h>
#endif 

#define DIETTAG 2

using namespace OpenSim;
using namespace SimTK;
using namespace std;

// Global variables
int evalCount = 0;
double bestSoFar = SimTK::Infinity;
double forceThreshold = 0.0;
double CLFThreshold = 0.0;

double maxActivation = 1.0;
double minActivation = -1.0;

double integratorTolerance = 1.0e-4;

// Optimizer settings and flags
int lambda = 100;
double sigma = 0.005;
int maxIterations = 1000;
bool doOpt = false;
bool doResume = false;

ofstream optLog("optLog.txt", ofstream::out);
ofstream fwdLog("fwdLog.txt", ofstream::out);
std::clock_t optimizerStartTime;

//=============================================================================
// Utility Methods
//=============================================================================
double evaluateModelAtFinalState(const Model& model, const State& s, bool outputLog, const Vector &newControls);

//=============================================================================
// OPTIMIZER SYSTEM: JumpingOptimizationSystem
// Defines a constructor and objective function for the optimization
//=============================================================================
class JumpingOptimizationSystem : public OptimizerSystem {
   public:

	   /* Constructor class. Parameters passed are accessed in the objectiveFunc() class. */
	   JumpingOptimizationSystem(int numParameters, State& s, Model& aModel, Manager& manager, Array<double> controllerTimePoints): 
             OptimizerSystem(numParameters),
			 si(s),
		     osimModel(aModel),
			 manager(manager),
			 nodeTimePoints(controllerTimePoints), 
			 numNodesPerActuator(controllerTimePoints.getSize())
	   {
		   // Get a copy of the state to reuse
		   _sCopy = s;
	   }
   
    void evalCtrl( const Vector &newControls, double &f ) const {
		std::clock_t startTime = std::clock(); 

        // reuse the initial state;
		_sCopy = si;
		State& s = _sCopy;


		// Grab actuator and controller sets
		int numControls = osimModel.getNumControls();
		PrescribedController* thisController
			= static_cast<PrescribedController*>(&osimModel.getControllerSet()[0]);

		// Update controller's functions
		FunctionSet& controlFunctions = thisController->upd_ControlFunctions();

		// Update the control values
		for (int i = 0; i < numControls; i++) {
			PiecewiseLinearFunction* pwlFunc
				= static_cast<PiecewiseLinearFunction*>(&controlFunctions[i]);
			if (pwlFunc){
				for (int j = 0; j < numNodesPerActuator; j++) {
					//cout << "line 150 j=" << j << endl;
					pwlFunc->setY(j, newControls[numNodesPerActuator*i + j]);
				}
			}

		}

		manager.integrate(s);
		f = evaluateModelAtFinalState(osimModel, s, false, newControls);
   }
			 	
	int objectiveFunc(const Vector &newControls, bool new_coefficients, Real& f ) const {
		evalCtrl(newControls, f);
		
		return 0;  
	}

private:
	Model& osimModel;
	Manager& manager;
	Array<double> nodeTimePoints;
	int numNodesPerActuator;
	// the initial state that all simulations must begin with
	const State& si;
	// keep a working copy of the state so each objective func eval does not have to make a copy
	mutable State _sCopy;
 };


//______________________________________________________________________________
/**
 * Define an optimization problem that finds a set of torque actuator controls to maximize 
 * the jumping distance with soft constraints.
 */
int main(int argc, char* argv[])
{
#ifdef USE_MPI
	int myrank; 
	int rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	int numtasks;
	numtasks = 0;
	myrank = 0;
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	printf ("Number of processes= %d My rank= %d\n", numtasks,myrank);
#endif
	std::clock_t startTime = std::clock();	

	string modelFile = "AshbyModel_twoConstraints.osim";
	string ctrlFile = "InitialControllerParameters.sto";

	// Command line options
	for (int i = 0; i < argc; i++) {
        
		if (strcmp(argv[i], "-m") == 0) {
			modelFile = argv[i+1];
		}

		if (strcmp(argv[i], "-cf") == 0) {
           	ctrlFile = argv[i+1];
		}

		if (strcmp(argv[i], "-resume") == 0) {
			doResume = true;
		}
        
		if (strcmp(argv[i], "-maxIters") == 0) {
			maxIterations = atoi(argv[i+1]);
		}

		if (strcmp(argv[i], "-lambda") == 0) {
			lambda = atoi(argv[i+1]);
		}

		if (strcmp(argv[i], "-sigma") == 0) {
			sigma = atof(argv[i+1]);
		}

		if (strcmp(argv[i], "-opt") == 0) {
			doOpt = true; 
		}
	}

	Model osimModel(modelFile);
	Storage initialControllerParametersFile = Storage(ctrlFile);
	Array<double> controllerNodeTimes; 
	initialControllerParametersFile.getTimeColumn(controllerNodeTimes);
	double initialTime = initialControllerParametersFile.getFirstTime();
	double finalTime = initialControllerParametersFile.getLastTime();

	int numActuators, numNodesPerActuator, numParameters;
	const Set<Actuator> &actuatorSet = osimModel.getActuators();

	// The number of parameters is (numActuators*numNodesPerActuator) for
    // piecewise linear functions of controls
	numActuators = actuatorSet.getSize();
	optLog << "numActuators = " << numActuators << endl;
	numNodesPerActuator = controllerNodeTimes.getSize();
	numParameters = numActuators*numNodesPerActuator;
	optLog << "numNodesPerActuator = " << numNodesPerActuator << endl;

	/* Define initial values for controllerParameters. Each controller is discretized by
	nodes as defined in InitialControlParameters.sto. */
		
	// initialize controller parameters
	Vector controllerParameters(numParameters, 0.05);
		
	PrescribedController* thisNewController = new PrescribedController();

	// Add prescribed controllers to each actuator and initialize controllerParameters.
	for (int i = 0; i < numActuators; i++) {
		Array<double> theseControllerParameters;
		initialControllerParametersFile.getDataColumn(i, theseControllerParameters);

		// Construct PiecewiseLinearFunction from the initial controller parameters
		PiecewiseLinearFunction thisPWLFunction;
		for (int j = 0; j < numNodesPerActuator; j++) {
			controllerParameters[numNodesPerActuator*i + j] = theseControllerParameters[j];
			//cout << "i=" << i << ", j=" << j << std::endl;
			thisPWLFunction.addPoint(controllerNodeTimes[j], controllerParameters[numNodesPerActuator*i + j]);
		}

		thisNewController->addActuator(actuatorSet.get(i));
		thisNewController->setName(actuatorSet.get(i).getName());
		thisNewController->prescribeControlForActuator(actuatorSet.get(i).getName(), thisPWLFunction.clone());
	}

	osimModel.addController(thisNewController);
	
    // Add necessary probes to track slipping and coordinate limit forces
	ConstraintSlipPenaltyProbe* toeConstraintSlipPenaltyProbe = new ConstraintSlipPenaltyProbe();
	toeConstraintSlipPenaltyProbe->setName("toeConstraintSlipPenaltyProbe");
	toeConstraintSlipPenaltyProbe->set_constraint_names(0, "toeConstraint");
	toeConstraintSlipPenaltyProbe->set_friction_coefficient(0.8);
	toeConstraintSlipPenaltyProbe->set_probe_operation("integrate");
	toeConstraintSlipPenaltyProbe->set_initial_conditions_for_integration(0, 0.0);
	osimModel.addProbe(toeConstraintSlipPenaltyProbe);

	ConstraintSlipPenaltyProbe* heelConstraintSlipPenaltyProbe = new ConstraintSlipPenaltyProbe();
	heelConstraintSlipPenaltyProbe->setName("heelConstraintSlipPenaltyProbe");
	heelConstraintSlipPenaltyProbe->set_constraint_names(0, "heelConstraint");
	heelConstraintSlipPenaltyProbe->set_friction_coefficient(0.8);
	heelConstraintSlipPenaltyProbe->set_probe_operation("integrate");
	heelConstraintSlipPenaltyProbe->set_initial_conditions_for_integration(0, 0.0);
	osimModel.addProbe(heelConstraintSlipPenaltyProbe);
		
	CoordinateLimitForceProbe* coordinateLimitForceProbe = new CoordinateLimitForceProbe();
	coordinateLimitForceProbe->setName("coordinateLimitForceProbe");
	coordinateLimitForceProbe->set_coordinate_limit_force_names(0, "all");
	coordinateLimitForceProbe->set_force_threshold(CLFThreshold);
	coordinateLimitForceProbe->set_exponent(2.0);
	coordinateLimitForceProbe->set_probe_operation("integrate");
	coordinateLimitForceProbe->set_initial_conditions_for_integration(0, 0.0);
	osimModel.addProbe(coordinateLimitForceProbe);

	// Initialize the system. initSystem() cannot be used here because adding the event handler
	// must be done between buildSystem() and initializeState().
	cout << "Building system and adding EventHandlers" << endl;
	osimModel.buildSystem();

    // Stop the simulation if the toe or heel hit the ground
	//TerminateSimulationToeHeight *terminateToeHeight = new TerminateSimulationToeHeight(osimModel, toeHeightThreshold);
	TerminateSimulationZeroVelAndAcc *terminateZeroVelAndAcc = new TerminateSimulationZeroVelAndAcc(osimModel, 1e-5, 1e-5);
	//TerminateSimulationHeelHeight *terminateHeelHeight = new TerminateSimulationHeelHeight(osimModel, heelHeightThreshold);

    // Stop the simulation if the vertical component of the toe or heel
    // constraint is pulling down on the model
	//TerminateSimulationToeForce *terminateToeForce = new TerminateSimulationToeForce(osimModel, 1.0);
	//TerminateSimulationHeelForce *terminateHeelForce = new TerminateSimulationHeelForce(osimModel, 1.0);

    // Detect when the vertical component of the toe or heel constraint reaches
    // 0 and release the constraint for flight
	//ReleaseToeConstraint *releaseToeConstraint = new ReleaseToeConstraint(osimModel, forceThreshold);
	ReleaseSeatConstraint *releaseSeatConstraint = new ReleaseSeatConstraint(osimModel, forceThreshold);
	//ReleaseHeelConstraint *releaseHeelConstraint = new ReleaseHeelConstraint(osimModel, forceThreshold);

    // Add event handlers to the model
	//osimModel.updMultibodySystem().addEventHandler(terminateToeHeight);
	//osimModel.updMultibodySystem().addEventHandler(terminateHeelHeight);
	//osimModel.updMultibodySystem().addEventHandler(terminateToeForce);
	//osimModel.updMultibodySystem().addEventHandler(terminateHeelForce);
	//osimModel.updMultibodySystem().addEventHandler(releaseToeConstraint);
	osimModel.updMultibodySystem().addEventHandler(releaseSeatConstraint);
	osimModel.updMultibodySystem().addEventHandler(terminateZeroVelAndAcc);
	//osimModel.updMultibodySystem().addEventHandler(releaseHeelConstraint);
		
	State &osimState = osimModel.initializeState();

	// Create the integrator for the simulation.
	RungeKuttaMersonIntegrator integrator(osimModel.getMultibodySystem());
	//SemiExplicitEuler2Integrator integrator(osimModel.getMultibodySystem());
	integrator.setAccuracy(integratorTolerance);

	// Create a manager to run the simulation. Can change manager options to save run time and memory or print more information
	Manager manager(osimModel, integrator);
	manager.setWriteToStorage(false);
	manager.setPerformAnalyses(false);

	// Integrate from initial time to final time and integrate
	manager.setInitialTime(initialTime);
	manager.setFinalTime(finalTime);	

	// Create the OptimizationSystem. Initialize the objective function value "f".
	JumpingOptimizationSystem sys(numParameters, osimState, osimModel, manager, controllerNodeTimes);
	Real f = NaN;
		
	// Set lower and upper bounds.
	Vector lower_bounds(numParameters, minActivation);
	Vector upper_bounds(numParameters, maxActivation);
		
	sys.setParameterLimits( lower_bounds, upper_bounds );
	optLog << "lower_bounds = " << lower_bounds << endl;
	optLog << "upper_bounds = " << upper_bounds << endl << endl;

	int failed = 0;
			
	if (doOpt) {
		CMAOptimizer opt(sys); optLog << "using CMA Optimizer" << endl;

		optLog << "maxIterations = " << maxIterations << endl;
		optLog << "lambda = " << lambda << endl;
		optLog << "sigma = " << sigma << endl;
		//optLog << "w1 = " << w1 << endl;
		//optLog << "w2 = " << w2 << endl;
		//optLog << "w3 = " << w3 << endl;
		//optLog << "w4 = " << w4 << endl;
		//optLog << "del_CGx = " << del_CGx << endl;
		//optLog << "del_CGy = " << del_CGy << endl;
		optLog << "integratorTolerance = " << integratorTolerance << endl;
		optLog << "mu = " << mu << endl;
		//optLog << "timePenaltyThreshold = " << timePenaltyThreshold << endl;
		optLog << "CLFThreshold = " << CLFThreshold << endl;

		opt.setMaxIterations(maxIterations);

		if (doResume) opt.setAdvancedBoolOption("resume", true);
		opt.setAdvancedIntOption("lambda", lambda);
		opt.setAdvancedRealOption("sigma", sigma);
#ifdef USE_MPI
		opt.setAdvancedBoolOption("enableMPI", true);
		if (numtasks == 1) opt.setAdvancedBoolOption("enableMPI", false);
#else
		opt.setAdvancedBoolOption("enableMPI", false);
#endif
		try{
			f = opt.optimize(controllerParameters);
		}
		catch (const std::exception& ex){
			std::cout << ex.what() << std::endl;
			failed = 1;
		}
	
		cout << "Elapsed time = " << (std::clock()-startTime)/CLOCKS_PER_SEC << "s" << endl;
		
		const Set<Actuator>& actuators = osimModel.getActuators();
		for(int i=0; i < actuators.getSize(); ++i){
			Vector thisControllerParameters(numNodesPerActuator, 0.0);
			for(int j=0; j < numNodesPerActuator; j++) {
				thisControllerParameters[j] = controllerParameters[i*numNodesPerActuator + j];
			}
		cout << actuators[i].getName() << " control values = " << thisControllerParameters << endl;
		}

		osimModel.print("optimizedModel.osim");

		std::string fileName = "optimizedControls.sto";
		std::ofstream optimizedControlsFile(fileName, std::ofstream::out);

		// HEADERS
		optimizedControlsFile << "optimizedControls" << endl;
		optimizedControlsFile << "version=1" << endl;
		optimizedControlsFile << "nRows=" << controllerParameters.size()/numActuators << endl;
		optimizedControlsFile << "nColumns=" << numActuators+1 << endl;
		optimizedControlsFile << "inDegrees=no" << endl;
		optimizedControlsFile << "endheader" << endl;

		optimizedControlsFile << "time" << "\t";
		optimizedControlsFile << "ankle_MLCA" << "\t";
		optimizedControlsFile << "knee_MLCA" << "\t";
		optimizedControlsFile << "hip_MLCA" << "\t";
		optimizedControlsFile << "shoulder_MLCA" << "\t";
		optimizedControlsFile << endl;
		optimizedControlsFile.precision(16);

		double timeStep = 0.05;
			
		for (int i = 0; i < controllerParameters.size()/numActuators; i++) {
			optimizedControlsFile << timeStep*i << "\t";
			for (int j = 0; j < numActuators; j++) {
				optimizedControlsFile << controllerParameters[i + j*controllerParameters.size()/numActuators] << "\t";
			}
			optimizedControlsFile << endl;
		}
		optimizedControlsFile.close();

		cout << "\nMinimum Objective Function Value = " << f << endl;

	}

    // Run a single forward simulation without optimization
    else {
        osimModel.print("fwdOnlyModel.osim");

        ForceReporter* forceReporter = new ForceReporter(&osimModel);
        forceReporter->setName("ForceReporter");
        forceReporter->includeConstraintForces(true);
        osimModel.addAnalysis(forceReporter);

        PointKinematics* pointKinematics_toe = new PointKinematics(&osimModel);
        pointKinematics_toe->setName("PointKinematics_TOE");
        pointKinematics_toe->setBodyPoint("foot", toeInFootPosition);
        osimModel.addAnalysis(pointKinematics_toe);

        PointKinematics* pointKinematics_heel = new PointKinematics(&osimModel);
        pointKinematics_heel->setName("PointKinematics_HEEL");
        pointKinematics_heel->setBodyPoint("foot", Vec3(0.0));
        osimModel.addAnalysis(pointKinematics_heel);

        startTime = std::clock();

        manager.setWriteToStorage(true);
        manager.setPerformAnalyses(true);

        // Integrate from initial time to final time and integrate
        manager.setInitialTime(initialTime);
        manager.setFinalTime(finalTime);
        manager.integrate(osimState);

        Storage states(manager.getStateStorage());
        states.print("states.sto");

        Storage forceStorage = forceReporter->getForceStorage();
        forceStorage.print("forceReporter.sto");
        pointKinematics_toe->getPositionStorage()->print("PointKinematics_TOE_pos.sto");
        pointKinematics_toe->getVelocityStorage()->print("PointKinematics_TOE_vel.sto");
        pointKinematics_toe->getAccelerationStorage()->print("PointKinematics_TOE_acc.sto");
        pointKinematics_heel->getPositionStorage()->print("PointKinematics_HEEL_pos.sto");
        pointKinematics_heel->getVelocityStorage()->print("PointKinematics_HEEL_vel.sto");
        pointKinematics_heel->getAccelerationStorage()->print("PointKinematics_HEEL_acc.sto");

        clock_t endTime = std::clock();
        fwdLog << "computeTime: " << endTime - startTime << "ms" << endl;

		evaluateModelAtFinalState(osimModel, osimState, true, controllerParameters);
	}

#ifdef USE_MPI
	for (int rank = 1; rank < numtasks; ++rank) {
		MPI_Send(0, 0, MPI_INT, rank, DIETTAG, MPI_COMM_WORLD);
	}
	MPI_Finalize();
#endif
	
	// End of main() routine.
	return failed;
}

double evaluateModelAtFinalState(const Model& osimModel, const State& s, bool outputLog, const Vector &newControls) {
	/* CALCULATING OBJECTIVE FUNCTION */
	osimModel.getMultibodySystem().realize(s, Stage::Position);

	// MAX DISTANCE
	Vec3 toeInGroundPosition, heelInGroundPosition, thighPosition, shankPosition, torsoPosition;
	osimModel.getSimbodyEngine().getPosition(s, osimModel.getBodySet().get("foot"), Vec3(0.0), heelInGroundPosition);
	osimModel.getSimbodyEngine().getPosition(s, osimModel.getBodySet().get("thigh"), Vec3(0.0), thighPosition);
	osimModel.getSimbodyEngine().getPosition(s, osimModel.getBodySet().get("shank"), Vec3(0.0), shankPosition);
	osimModel.getSimbodyEngine().getPosition(s, osimModel.getBodySet().get("torso"), Vec3(0.0), torsoPosition);
	Vec3 torsoVelocityInOrigin;
	osimModel.getSimbodyEngine().getVelocity(s, osimModel.getBodySet().get("torso"), origin, torsoVelocityInOrigin);
	Vec3 torsoAccelerationInOrigin;
	osimModel.getSimbodyEngine().getAcceleration(s, osimModel.getBodySet().get("torso"), origin, torsoAccelerationInOrigin);

	//Zero velocity and acceleration at final state
	double velocityPenalty = torsoVelocityInOrigin.scalarNormSqr();
	double accelerationPenalty = torsoAccelerationInOrigin.scalarNormSqr();
	
	
	//double maxDistance = heelInGroundPosition[0];
	double K_CGy = 0.0;
	double K_CGx = 0.0;
	Vec3 COM_position = osimModel.calcMassCenterPosition(s);
	double ry = COM_position[1];
	double rx = COM_position[0];

	double del_CGx = 0.0; // horizontal position landing constraint
	double del_CGy = 0.9; // vertical position landing constraint


	// STANDING CONSTRAINTS
	if (ry - heelInGroundPosition[1] < del_CGy) {
		K_CGy = heelInGroundPosition[1] - ry + del_CGy; //if less than del_CGy, add to f (makes it try to reach del_CGy)
	}
	K_CGx = abs(heelInGroundPosition[0] - rx - del_CGx);
	double alignmentError = pow(thighPosition[0] - heelInGroundPosition[0], 2.0) + pow(shankPosition[0] - heelInGroundPosition[0], 2.0) + pow(torsoPosition[0] - heelInGroundPosition[0], 2.0);

	//double standingError = w1*K_CGy + alignmentError;

	const ProbeSet probeSet = osimModel.getProbeSet();

	// SLIP CONSTRAINT
	const ConstraintSlipPenaltyProbe* toeSlipProbe = dynamic_cast<const ConstraintSlipPenaltyProbe*>(&probeSet.get("toeConstraintSlipPenaltyProbe"));
	SimTK::Vector toeSlipProbeOutput = toeSlipProbe->getProbeOutputs(s);
	const ConstraintSlipPenaltyProbe* heelSlipProbe = dynamic_cast<const ConstraintSlipPenaltyProbe*>(&probeSet.get("heelConstraintSlipPenaltyProbe"));
	SimTK::Vector heelSlipProbeOutput = heelSlipProbe->getProbeOutputs(s);

	// COORDINATE LIMIT FORCE CONSTRAINT //if use ligaments, penalizes you.
	const CoordinateLimitForceProbe* coordinateLimitForceProbe = dynamic_cast<const CoordinateLimitForceProbe*>(&probeSet.get("coordinateLimitForceProbe"));
	SimTK::Vector coordinateLimitForceProbeOutput = coordinateLimitForceProbe->getProbeOutputs(s);

	// TIME PENALTY
	double timePenalty = 0.0;
	timePenalty = s.getTime();

	// time penalty
	double timePenaltyThreshold = 1.2;
	//if (s.getTime() < timePenaltyThreshold) {
	//	timePenalty += 1 / (s.getTime()*s.getTime()) - 1 / (timePenaltyThreshold*timePenaltyThreshold);
	//}


	// Actuator Controls Usage Penalty (optimize for less output forces)
	double MLCAControlsUsage = 0.0;
	for (int i = 0; i < newControls.size(); i++) {
		MLCAControlsUsage = MLCAControlsUsage + pow(newControls(i), 2.0);
	}

	// variables related to objective function
	double w0 = 1.0; // alignment error constraint weight
	double w1 = 10.0; // height constraint weight
	double w2 = 0.0; //1.0e-4; // ligament penalty weight
	double w3 = 0.1; // MLCA actuator weight
	double w4 = 0.0; // 10.0; // time penalty weight


	double f = w0*alignmentError + w1*K_CGy + w2*coordinateLimitForceProbeOutput(0) + w3*MLCAControlsUsage + w4*timePenalty;
	//w1*(K_CGy + K_CGx); // + w2*coordinateLimitForceProbeOutput(0) + w3*(toeSlipProbeOutput(0) + heelSlipProbeOutput(0)) + w4*timePenalty;


	if (outputLog) {
		fwdLog << "endOfSimulationTime: " << s.getTime() << endl;
		fwdLog << "integratorTolerance: " << integratorTolerance << endl;

		fwdLog << "f = " << f << endl;

		fwdLog << "w0 = " << w0 << ", alignmentError = " << alignmentError << ", contribution = " << w0*alignmentError / f << endl;

		fwdLog << "w1 = " << w1 << ", K_CGy = " << K_CGy << ", contribution = " << w1*K_CGy / f << endl;

		fwdLog << "w2 = " << w2 << ", coordinateLimitForceProbeOutput = " << coordinateLimitForceProbeOutput << ", contribution = " << w2*coordinateLimitForceProbeOutput(0) / f << endl;

		fwdLog << "w3 = " << w3 << ", MLCA Control Usage = " << MLCAControlsUsage << ", contribution = " << w3*MLCAControlsUsage / f << endl;

		fwdLog << "w4 = " << w4 << ", timePenalty = " << timePenalty << ", contribution = " << w4*timePenalty / f << endl;

		//fwdLog << "toeSlipProbeOutput" << toeSlipProbeOutput << endl;
		//fwdLog << "heelSlipProbeOutput" << heelSlipProbeOutput << endl;
	}

	static double bestf = 1e10;

	if (f < bestf){
		bestf = f;
		cout << "f = " << f;
		cout << ", align=" << w0*alignmentError;
		cout << ", height=" << w1*K_CGy;
		cout << ", injury=" << w2*coordinateLimitForceProbeOutput(0);
		cout << ", MLCA= " << w3*MLCAControlsUsage;
		cout << ", time=" << w4*timePenalty;
		cout << ", endTime=" << s.getTime();
		cout << ", vel=" << velocityPenalty;
		cout << ", acc=" << accelerationPenalty << endl;
	}

	return f;

}