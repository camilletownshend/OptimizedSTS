#include <OpenSim/OpenSim.h>

using namespace OpenSim;
using namespace SimTK;
using namespace std;

double toeHeightThreshold = 0.0;
double heelHeightThreshold = 0.0;
double seatHeightThreshold = 0.0;
Vec3 toeInFootPosition = Vec3(0.14, 0.0, 0.0);
Vec3 origin = Vec3(0.0, 0.0, 0.0);
double mu = 0.8;
//double timeThreshold = 1.2;

double toeConstraintTimeOff = 0.0;
double heelConstraintTimeOff = 0.0;
double seatConstraintTimeOff = 0.0;

double GRFLowestLocalMinimum = 1.0;
double systemMass = 709;

//=============================================================================
// EVENT HANDLER: TerminateSimulationZeroVelAndAcc
// Triggered when the model is stationary.
//=============================================================================
class TerminateSimulationZeroVelAndAcc : public SimTK::TriggeredEventHandler {
public:

	// CONSTRUCTOR
	TerminateSimulationZeroVelAndAcc(const Model& m, double velThreshold, double accThreshold) : TriggeredEventHandler(Stage::Acceleration),
		_model(m), _velThreshold(velThreshold), _accThreshold(accThreshold)
	{
		getTriggerInfo().setTriggerOnRisingSignTransition(false);
		getTriggerInfo().setTriggerOnFallingSignTransition(true);
	}

	// WITNESS FUNCTION
	SimTK::Real getValue(const SimTK::State& s) const
	{
		
		Vec3 torsoVelocityInOrigin;
		_model.getSimbodyEngine().getVelocity(s, _model.getBodySet().get("torso"), origin, torsoVelocityInOrigin);
		
		Vec3 torsoAccelerationInOrigin;
		_model.getSimbodyEngine().getAcceleration(s, _model.getBodySet().get("torso"), origin, torsoAccelerationInOrigin);

		//cout << "time = " << s.getTime() << ", toeEvent = " << torsoVelocityInOrigin[1] - toeHeightThreshold + 1.0e-5 << endl;

		if (torsoVelocityInOrigin.scalarNormSqr() < _velThreshold && torsoAccelerationInOrigin.scalarNormSqr() < _accThreshold)
			return 0.0;
		
		return 1;
	}

	// EVENT HANDLER FUNCTION
	void handleEvent(SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
	{
		//GRFLowestLocalMinimum = 1.0;
		//if (s.getTime() > timeThreshold) {
		terminate = true;
		//}
		//cout << "simulation terminated on toe height" << endl;
		//cout << "toeHeightThreshold: " << toeHeightThreshold << endl;

	}


private:
	const Model& _model;
	double _velThreshold;
	double _accThreshold;
};



//=============================================================================
// EVENT HANDLER: TerminateSimulationToeHeight
// Triggered when the GRF normal falls below some threshold
//=============================================================================
class TerminateSimulationToeHeight : public SimTK::TriggeredEventHandler {
public:
    
    // CONSTRUCTOR
    TerminateSimulationToeHeight(const Model& m, double threshold) : TriggeredEventHandler(Stage::Acceleration), 
        _model(m), _heightThreshold(threshold)
    {
		getTriggerInfo().setTriggerOnRisingSignTransition(false);
		getTriggerInfo().setTriggerOnFallingSignTransition(true);
    }

    // WITNESS FUNCTION
    SimTK::Real getValue(const SimTK::State& s) const
	{
		// keep returning a positive value while constraint is active
		if (_model.getConstraintSet().getSize() > 0){
			if (!_model.getConstraintSet().get("toeConstraint").isDisabled(s)) {
				return 1.0;
			}
		}

		// minimum of toe and heel position
		Vec3 toeInGroundPosition;
		_model.getSimbodyEngine().getPosition(s, _model.getBodySet().get("foot"), toeInFootPosition, toeInGroundPosition);

		//cout << "time = " << s.getTime() << ", toeEvent = " << toeInGroundPosition[1] - toeHeightThreshold + 1.0e-5 << endl;

		return toeInGroundPosition[1] - toeHeightThreshold + 3.0e-2;

	}

    // EVENT HANDLER FUNCTION
    void handleEvent (SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
	{
		//GRFLowestLocalMinimum = 1.0;
		//if (s.getTime() > timeThreshold) {
			terminate = true;
		//}
		//cout << "simulation terminated on toe height" << endl;
		//cout << "toeHeightThreshold: " << toeHeightThreshold << endl;

	}


private:
    const Model& _model;
    double _heightThreshold;
};

//=============================================================================
// EVENT HANDLER: TerminateSimulationHeelHeight
// Triggered when the GRF normal falls below some threshold
//=============================================================================
class TerminateSimulationHeelHeight : public SimTK::TriggeredEventHandler {
public:
    
    // CONSTRUCTOR
    TerminateSimulationHeelHeight(const Model& m, double threshold) : TriggeredEventHandler(Stage::Acceleration), 
        _model(m), _heightThreshold(threshold)
    {
		getTriggerInfo().setTriggerOnRisingSignTransition(false);
		getTriggerInfo().setTriggerOnFallingSignTransition(true);
    }

    // WITNESS FUNCTION
    SimTK::Real getValue(const SimTK::State& s) const
	{
		if (_model.getConstraintSet().getSize() > 0){
			// keep returning a positive value while constraint is active
			if (!_model.getConstraintSet().get("heelConstraint").isDisabled(s)) {
				return 1.0;
			}
		}

		// minimum of toe and heel position
		Vec3 heelInGroundPosition;
		_model.getSimbodyEngine().getPosition(s, _model.getBodySet().get("foot"), Vec3(0.0), heelInGroundPosition);

		return heelInGroundPosition[1] - heelHeightThreshold + 3.0e-2;

	}

    // EVENT HANDLER FUNCTION
    void handleEvent (SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
	{
		//GRFLowestLocalMinimum = 1.0;
		//if (s.getTime() > timeThreshold) {
			terminate = true;
		//}
		//cout << "simulation terminated on heel height" << endl;
		//cout << "heelHeightThreshold: " << heelHeightThreshold << endl;

	}


private:
    const Model& _model;
    double _heightThreshold;
};

//=============================================================================
// EVENT HANDLER: TeminateSimulationToeForce
// Triggered when the GRF normal falls below some threshold
//=============================================================================
class TerminateSimulationToeForce : public SimTK::TriggeredEventHandler {
public:
    
    // CONSTRUCTOR
	TerminateSimulationToeForce(Model& m, double forceThreshold = 0.0) : TriggeredEventHandler(Stage::Acceleration),
		_model(m), _forceThreshold(forceThreshold)
    {
		getTriggerInfo().setTriggerOnRisingSignTransition(false);
		getTriggerInfo().setTriggerOnFallingSignTransition(true);
    }

    // WITNESS FUNCTION
    SimTK::Real getValue(const SimTK::State& s) const
	{	
		SimTK::Real verticalForceOnFootFromGround = 0;

		if (_model.getConstraintSet().getSize() > 0){
			const OpenSim::Constraint &toeConstraint = _model.getConstraintSet().get("toeConstraint");

			if (toeConstraint.isDisabled(s)) {
				return -1.0;
			}

			SimTK::ConstraintIndex toeConstraintIndex(0);
			SimTK::Constraint& simConstraint = _model.updMatterSubsystem().updConstraint(toeConstraintIndex);

			// number of bodies being directly constrained
			int ncb = simConstraint.getNumConstrainedBodies();
			// number of mobilities being directly constrained
			int ncm = simConstraint.getNumConstrainedU(s);

			SimTK::Vector_<SimTK::SpatialVec> bodyForcesInAncestor(ncb);
			bodyForcesInAncestor.setToZero();
			SimTK::Vector mobilityForces(ncm, 0.0);

			toeConstraint.calcConstraintForces(s, bodyForcesInAncestor, mobilityForces);
			verticalForceOnFootFromGround = bodyForcesInAncestor(0)[1][1] + _forceThreshold;
		}
		else{
			//const OpenSim::Force& floorContact = _model.getForceSet().get("foot_floor");
			const OpenSim::Force& floorContact_r = _model.getForceSet().get("foot_r");
			const OpenSim::Force& floorContact_l = _model.getForceSet().get("foot_l");
			//Array<double> grfs = floorContact.getRecordValues(s);
			Array<double> grfs_r = floorContact_r.getRecordValues(s);
			Array<double> grfs_l = floorContact_l.getRecordValues(s);
			//verticalForceOnFootFromGround = grfs[1]+_forceThreshold; // if you want to include flight
			verticalForceOnFootFromGround = -grfs_r[1] - grfs_l[1]; // take-off only
		}


		return verticalForceOnFootFromGround;

	}

    // EVENT HANDLER FUNCTION
    void handleEvent (SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
	{
		//GRFLowestLocalMinimum = 1.0;
		terminate = true;
		//cout << "Terminate simulation on toe force" << endl;
	}


private:
    Model& _model;
	double _forceThreshold;
};

//=============================================================================
// EVENT HANDLER: TeminateSimulationHeelForce
// Triggered when the GRF normal falls below some threshold
//=============================================================================
class TerminateSimulationHeelForce : public SimTK::TriggeredEventHandler {
public:
    
    // CONSTRUCTOR
    TerminateSimulationHeelForce(Model& model, double forceThreshold=0.0) : 
		TriggeredEventHandler(Stage::Acceleration), 
		_model(model), _forceThreshold(forceThreshold)
    {
		getTriggerInfo().setTriggerOnRisingSignTransition(true);
		getTriggerInfo().setTriggerOnFallingSignTransition(false);
    }

    // WITNESS FUNCTION
    SimTK::Real getValue(const SimTK::State& s) const
	{	
		SimTK::Real verticalForceOnFootFromGround = 0;

		if (_model.getConstraintSet().getSize() > 0){
			const OpenSim::Constraint &heelConstraint = _model.getConstraintSet().get("heelConstraint");

			if (heelConstraint.isDisabled(s)) {
				return -1.0;
			}

			SimTK::ConstraintIndex heelConstraintIndex(1);
			SimTK::Constraint& simConstraint = _model.updMatterSubsystem().updConstraint(heelConstraintIndex);

			// number of bodies being directly constrained
			int ncb = simConstraint.getNumConstrainedBodies();
			// number of mobilities being directly constrained
			int ncm = simConstraint.getNumConstrainedU(s);

			SimTK::Vector_<SimTK::SpatialVec> bodyForcesInAncestor(ncb);
			bodyForcesInAncestor.setToZero();
			SimTK::Vector mobilityForces(ncm, 0.0);

			heelConstraint.calcConstraintForces(s, bodyForcesInAncestor, mobilityForces);
			verticalForceOnFootFromGround = bodyForcesInAncestor(0)[1][1] + _forceThreshold;
		}
		else{
			const OpenSim::Force& floorContact = _model.getForceSet().get("foot_floor");
			Array<double> grfs = floorContact.getRecordValues(s);
			verticalForceOnFootFromGround = -grfs[1];
		}
		
		return verticalForceOnFootFromGround-_forceThreshold;
	}

    // EVENT HANDLER FUNCTION
    void handleEvent (SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
	{
		//GRFLowestLocalMinimum = 1.0;
		terminate = true;
		//cout << "Terminate simulation on heel force" << endl;
	}


private:
    Model& _model;
	double _forceThreshold;
};

//=============================================================================
// EVENT HANDLER: ReleaseToeConstraint
// Triggered when the GRF normal falls below some threshold
//=============================================================================
class ReleaseToeConstraint : public SimTK::TriggeredEventHandler {
public:
    
    // CONSTRUCTOR
    ReleaseToeConstraint(Model& m, double threshold) : TriggeredEventHandler(Stage::Acceleration), 
        _model(m), _forceThreshold(threshold)
    {
		getTriggerInfo().setTriggerOnRisingSignTransition(false);
		getTriggerInfo().setTriggerOnFallingSignTransition(true);
    }

    // WITNESS FUNCTION
    SimTK::Real getValue(const SimTK::State& s) const
	{

		const OpenSim::Constraint &toeConstraint = _model.getConstraintSet().get("toeConstraint");

		if (toeConstraint.isDisabled(s)) {
			return -1.0;
		}
		
		SimTK::ConstraintIndex toeConstraintIndex(0);
		SimTK::Constraint& simConstraint = _model.updMatterSubsystem().updConstraint(toeConstraintIndex);

		// number of bodies being directly constrained
		int ncb = simConstraint.getNumConstrainedBodies();
		// number of mobilities being directly constrained
		int ncm = simConstraint.getNumConstrainedU(s);

		SimTK::Vector_<SimTK::SpatialVec> bodyForcesInAncestor(ncb);
		bodyForcesInAncestor.setToZero();
		SimTK::Vector mobilityForces(ncm, 0.0);

		toeConstraint.calcConstraintForces(s, bodyForcesInAncestor, mobilityForces);

		//cout << "bodyForcesInAncestor: " << bodyForcesInAncestor << endl;
		//cout << "mobilityForces: " << mobilityForces << endl;

		Array<double> values(0.0,6*ncb+ncm);
		for(int i=0; i<ncb; ++i){
			for(int j=0; j<3; ++j){
				// Simbody constraints have reaction moments first and OpenSim reports forces first
				// so swap them here
				values[i*6+j] = (bodyForcesInAncestor(i)[1])[j]; // moments on constrained body i
				values[i*6+3+j] = (bodyForcesInAncestor(i)[0])[j]; // forces on constrained body i
			}
		}
		
		for(int i=0; i<ncm; ++i){
			values[6*ncb+i] = mobilityForces[i];
		}

		//cout << toeConstraint.getRecordLabels() << endl;
		//cout << values << endl;
		//cout << endl;
		double forceOnFootFromGround = values[1];
		//cout << values[1] << endl;

		return forceOnFootFromGround - _forceThreshold;

	}

    // EVENT HANDLER FUNCTION
    void handleEvent (SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
	{
		_model.updConstraintSet().get("toeConstraint").setDisabled(s, true);

		_model.getMultibodySystem().realize(s, SimTK::Stage::Position);
		Vec3 toeInGroundPosition;
		_model.getSimbodyEngine().getPosition(s, _model.getBodySet().get("foot"), toeInFootPosition, toeInGroundPosition);

		toeHeightThreshold = toeInGroundPosition[1];

		toeConstraintTimeOff = s.getTime();

		//cout << "toe constraint removed due to GRF going negative" << endl; //debug line
		//cout << "toeHeightThreshold set to " << toeHeightThreshold << endl;
		//system("pause");

		//cout << "toe constraint disabled at time = " << toeConstraintTimeOff << endl;

		//Vec3 COM_pos = _model.getMatterSubsystem().calcSystemMassCenterLocationInGround(s);
		//cout << "COM_pos = " << COM_pos << endl;
		//Vec3 COM_vel = _model.getMatterSubsystem().calcSystemMassCenterVelocityInGround(s);
		//cout << "COM_vel" << COM_vel << endl;

		

		terminate = false;
	}


private:
    Model& _model;
    double _forceThreshold;
};

//=============================================================================
// EVENT HANDLER: ReleaseHeelConstraint
// Triggered when the GRF normal falls below some threshold
//=============================================================================
class ReleaseHeelConstraint : public SimTK::TriggeredEventHandler {
public:
    
    // CONSTRUCTOR
    ReleaseHeelConstraint(Model& m, double threshold) : TriggeredEventHandler(Stage::Acceleration), 
        _model(m), _forceThreshold(threshold)
    {
		getTriggerInfo().setTriggerOnRisingSignTransition(false);
		getTriggerInfo().setTriggerOnFallingSignTransition(true);
    }

    // WITNESS FUNCTION
    SimTK::Real getValue(const SimTK::State& s) const
	{	

		const OpenSim::Constraint &heelConstraint = _model.getConstraintSet().get("heelConstraint");

		if (heelConstraint.isDisabled(s)) {
			return -1.0;
		}
		
		SimTK::ConstraintIndex heelConstraintIndex(1);
		SimTK::Constraint& simConstraint = _model.updMatterSubsystem().updConstraint(heelConstraintIndex);

		// number of bodies being directly constrained
		int ncb = simConstraint.getNumConstrainedBodies();
		// number of mobilities being directly constrained
		int ncm = simConstraint.getNumConstrainedU(s);

		SimTK::Vector_<SimTK::SpatialVec> bodyForcesInAncestor(ncb);
		bodyForcesInAncestor.setToZero();
		SimTK::Vector mobilityForces(ncm, 0.0);

		heelConstraint.calcConstraintForces(s, bodyForcesInAncestor, mobilityForces);

		Array<double> values(0.0,6*ncb+ncm);
		for(int i=0; i<ncb; ++i){
			for(int j=0; j<3; ++j){
				// Simbody constraints have reaction moments first and OpenSim reports forces first
				// so swap them here
				values[i*6+j] = (bodyForcesInAncestor(i)[1])[j]; // moments on constrained body i
				values[i*6+3+j] = (bodyForcesInAncestor(i)[0])[j]; // forces on constrained body i
			}
		}
		
		for(int i=0; i<ncm; ++i){
			values[6*ncb+i] = mobilityForces[i];
		}

		double forceOnHeelFromGround = values[1];

		return forceOnHeelFromGround - _forceThreshold;

	}

    // EVENT HANDLER FUNCTION
    void handleEvent (SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
	{
		_model.updConstraintSet().get("heelConstraint").setDisabled(s, true);

		_model.getMultibodySystem().realize(s, SimTK::Stage::Position);
		Vec3 heelInGroundPosition;
		_model.getSimbodyEngine().getPosition(s, _model.getBodySet().get("foot"), Vec3(0.0, 0.0, 0.0), heelInGroundPosition);

		heelHeightThreshold = heelInGroundPosition[1];

		heelConstraintTimeOff = s.getTime();

		//cout << "heel constraint removed due to GRF going negative" << endl; //debug line
		//cout << "heelHeightThreshold set to " << heelHeightThreshold << endl;
		//system("pause");
		terminate = false;
	}


private:
    Model& _model;
    double _forceThreshold;
};

//=============================================================================
// EVENT HANDLER: ReleaseSeatConstraint
// Triggered when the GRF normal falls below some threshold
//=============================================================================
class ReleaseSeatConstraint : public SimTK::TriggeredEventHandler {
public:

	// CONSTRUCTOR
	ReleaseSeatConstraint(Model& m, double threshold) : TriggeredEventHandler(Stage::Acceleration),
		_model(m), _forceThreshold(threshold)
	{
		getTriggerInfo().setTriggerOnRisingSignTransition(false);
		getTriggerInfo().setTriggerOnFallingSignTransition(true);
	}

	// WITNESS FUNCTION
	SimTK::Real getValue(const SimTK::State& s) const
	{

		const OpenSim::Constraint &seatConstraint = _model.getConstraintSet().get("seatConstraint");

		if (seatConstraint.isDisabled(s)) {
			return -1.0;
		}

		SimTK::ConstraintIndex seatConstraintIndex(1);
		SimTK::Constraint& simConstraint = _model.updMatterSubsystem().updConstraint(seatConstraintIndex);

		// number of bodies being directly constrained
		int ncb = simConstraint.getNumConstrainedBodies();
		// number of mobilities being directly constrained
		int ncm = simConstraint.getNumConstrainedU(s);

		SimTK::Vector_<SimTK::SpatialVec> bodyForcesInAncestor(ncb);
		bodyForcesInAncestor.setToZero();
		SimTK::Vector mobilityForces(ncm, 0.0);

		seatConstraint.calcConstraintForces(s, bodyForcesInAncestor, mobilityForces);

		Array<double> values(0.0, 6 * ncb + ncm);
		for (int i = 0; i<ncb; ++i){
			for (int j = 0; j<3; ++j){
				// Simbody constraints have reaction moments first and OpenSim reports forces first
				// so swap them here
				values[i * 6 + j] = (bodyForcesInAncestor(i)[1])[j]; // moments on constrained body i
				values[i * 6 + 3 + j] = (bodyForcesInAncestor(i)[0])[j]; // forces on constrained body i
			}
		}

		for (int i = 0; i<ncm; ++i){
			values[6 * ncb + i] = mobilityForces[i];
		}

		//std::cout << s.getTime() << " " << values << std::endl;
		double forceOnSeatFromGround = values[1];
		//cout << "forceOnSeatFromGround= " << values[1] << endl;
		return forceOnSeatFromGround - _forceThreshold;

	}

	// EVENT HANDLER FUNCTION
	void handleEvent(SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
	{
		_model.updConstraintSet().get("seatConstraint").setDisabled(s, true);

		_model.getMultibodySystem().realize(s, SimTK::Stage::Position);
		Vec3 seatInGroundPosition;
		_model.getSimbodyEngine().getPosition(s, _model.getBodySet().get("thigh"), Vec3(0.0, 0.0, 0.0), seatInGroundPosition);

		seatHeightThreshold = seatInGroundPosition[1];

		seatConstraintTimeOff = s.getTime();

		//cout << "seat constraint removed due to GRF going negative" << endl; //debug line
		//cout << "seatHeightThreshold set to " << seatHeightThreshold << endl;
		//system("pause");
		terminate = false;
	}


private:
	Model& _model;
	double _forceThreshold;
};

//=============================================================================
// EVENT HANDLER: GRFLowestLocalMinimumPenalty
// 
//=============================================================================
class GRFLowestLocalMinimumPenalty : public SimTK::TriggeredEventHandler {
public:
    
    // CONSTRUCTOR
    GRFLowestLocalMinimumPenalty(Model& m, double velocityThreshold) : TriggeredEventHandler(Stage::Acceleration), 
        _model(m), _velocityThreshold(velocityThreshold)
    {
		getTriggerInfo().setTriggerOnRisingSignTransition(true);
		getTriggerInfo().setTriggerOnFallingSignTransition(false);
    }

    // WITNESS FUNCTION
    SimTK::Real getValue(const SimTK::State& s) const
	{	

		Vec3 COM_Velocity = _model.getMatterSubsystem().calcSystemMassCenterVelocityInGround(s);
		
		return COM_Velocity[1] - _velocityThreshold;
		
	}

    // EVENT HANDLER FUNCTION
    void handleEvent (SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
	{
		//cout << "checking minimum GRF so far " << endl;
		SimTK::Vector verticalGRFProbeOutput = _model.getProbeSet().get("verticalGRFProbe").getProbeOutputs(s);
		
		GRFLowestLocalMinimum = verticalGRFProbeOutput(0);

		//cout << "GRFLowestLocalMinimum is now = " << GRFLowestLocalMinimum << " at time = " << s.getTime() << endl;

		terminate = false;
	}


private:
    Model& _model;
    double _velocityThreshold;
};


//=============================================================================
// EVENT REPORTER: SlipPenaltyReporter
// Keep track of penalty for slipping. Must be passed to OptimizationSystem to
// access _slipPenalty.
//=============================================================================
class SlipPenaltyReporter : public SimTK::PeriodicEventReporter {
public:
	SlipPenaltyReporter(Model& m,  Real interval) 
    :   PeriodicEventReporter(interval), _model(m), _interval(interval) {
		_slipPenalty = 0.0;
	}

    void handleEvent(const State& s) const {
	    const OpenSim::Constraint &toeConstraint = _model.getConstraintSet().get("toeConstraint");
		const OpenSim::Constraint &heelConstraint = _model.getConstraintSet().get("heelConstraint");

		SimTK::ConstraintIndex toeConstraintIndex(0);
		SimTK::Constraint& simConstraint = _model.updMatterSubsystem().updConstraint(toeConstraintIndex);

		// number of bodies being directly constrained
		int ncb = simConstraint.getNumConstrainedBodies();
		// number of mobilities being directly constrained
		int ncm = simConstraint.getNumConstrainedU(s);

		SimTK::Vector_<SimTK::SpatialVec> bodyForcesInAncestor(ncb);
		bodyForcesInAncestor.setToZero();
		SimTK::Vector mobilityForces(ncm, 0.0);

		Array<double> toeValues(0.0,6*ncb+ncm);
		Array<double> heelValues(0.0,6*ncb+ncm);

		if (!toeConstraint.isDisabled(s)) {
		
			toeConstraint.calcConstraintForces(s, bodyForcesInAncestor, mobilityForces);

			for(int i=0; i<ncb; ++i){
				for(int j=0; j<3; ++j){
					// Simbody constraints have reaction moments first and OpenSim reports forces first
					// so swap them here
					toeValues[i*6+j] = (bodyForcesInAncestor(i)[1])[j]; // moments on constrained body i
					toeValues[i*6+3+j] = (bodyForcesInAncestor(i)[0])[j]; // forces on constrained body i
				}
			}
		
			for(int i=0; i<ncm; ++i){
				toeValues[6*ncb+i] = mobilityForces[i];
			}
		}
		
		if (!heelConstraint.isDisabled(s)) {
		
			heelConstraint.calcConstraintForces(s, bodyForcesInAncestor, mobilityForces);

			for(int i=0; i<ncb; ++i){
				for(int j=0; j<3; ++j){
					// Simbody constraints have reaction moments first and OpenSim reports forces first
					// so swap them here
					heelValues[i*6+j] = (bodyForcesInAncestor(i)[1])[j]; // moments on constrained body i
					heelValues[i*6+3+j] = (bodyForcesInAncestor(i)[0])[j]; // forces on constrained body i
				}
			}
		
			for(int i=0; i<ncm; ++i){
				heelValues[6*ncb+i] = mobilityForces[i];
			}


		}


		if (abs(toeValues[0]) > mu*toeValues[1]) {
			_slipPenalty += (abs(toeValues[0]) - mu*toeValues[1])*_interval;
		}

		if (abs(heelValues[0]) > mu*heelValues[1]) {
			_slipPenalty += (abs(heelValues[0]) - mu*heelValues[1])*_interval;
		}

    }

	double getSlipPenalty() {
		return _slipPenalty;
	}

	
private:
    Model& _model;
	mutable double _slipPenalty; // seriously doesn't change the model or state
	double _interval;
};