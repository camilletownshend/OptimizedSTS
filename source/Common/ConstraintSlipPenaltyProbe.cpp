/* -------------------------------------------------------------------------- *
 *                   OpenSim:  ConstraintSlipPenaltyProbe.cpp                 *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2014 Stanford University and the Authors                *
 * Author(s): Carmichael Ong                                                  *
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

//=============================================================================
// INCLUDES
//=============================================================================
#include "ConstraintSlipPenaltyProbe.h"
#include <OpenSim/Simulation/Model/ConstraintSet.h>

using namespace std;
using namespace SimTK;
using namespace OpenSim;


//=============================================================================
// CONSTRUCTOR(S) AND SETUP
//=============================================================================
//_____________________________________________________________________________
/**
 * Default constructor.
 */
ConstraintSlipPenaltyProbe::ConstraintSlipPenaltyProbe() 
{
    setNull();
    constructProperties();
}


//_____________________________________________________________________________
/**
 * Set the data members of this TaskTermProbe to their null values.
 */
void ConstraintSlipPenaltyProbe::setNull(void)
{
    setAuthors("Carmichael Ong");
    _constraintIndex.clear();
}

//_____________________________________________________________________________
/**
 * Connect properties to local pointers.
 */
void ConstraintSlipPenaltyProbe::constructProperties(void)
{
    constructProperty_constraint_names();
    constructProperty_friction_coefficient(0.8);
}



//=============================================================================
// MODEL COMPONENT METHODS
//=============================================================================
//_____________________________________________________________________________
/**
 * Perform some set up functions that happen after the
 * object has been deserialized or copied.
 *
 * @param model OpenSim model containing this MuscleActivationProbe.
 */
void ConstraintSlipPenaltyProbe::connectToModel(Model& model)
{
    Super::connectToModel(model);
    _constraintIndex.clear();

    // Check to see if 'all' constraints are selected for probing.
    if (get_constraint_names(0) == "all") {
        Array<string> allConstraintNames;
		_model->getConstraintSet().getNames(allConstraintNames);
        set_constraint_names(allConstraintNames);
        //cout << "Set to all constraints: " << allConstraintNames << endl;
    }

    // Check that each Constraint in the constraint_names array exists in the model.
    // If so, populate the _constraintIndex internal array.
    const int numConstraints = getProperty_constraint_names().size();
    for (int i=0; i<numConstraints; i++) {
        const string& constraintName = get_constraint_names(i);
        const int k = model.getConstraintSet().getIndex(constraintName);
        if (k<0) {
            string errorMessage = getConcreteClassName() + ": Invalid Constraint '" 
                + constraintName + "' specified in <constraint_names>.";
            throw (Exception(errorMessage.c_str()));
        }
        else
            _constraintIndex.push_back(k);
    }

    // Sanity check.
    if (numConstraints != int(_constraintIndex.size()))
        throw (Exception("Size of _constraintIndex does not match number of constraints."));
}





//=============================================================================
// COMPUTATION
//=============================================================================
//_____________________________________________________________________________
/**
 * Compute the summed muscle activations.
 */
Vector ConstraintSlipPenaltyProbe::computeProbeInputs(const State& s) const
{
	//cout << "time: " << s.getTime() << endl;
	
	const int numConstraints = getProperty_constraint_names().size();
    const Set<Constraint>& cs = _model->getConstraintSet();
    SimTK::Vector TotalPenalty(1);
    TotalPenalty = 0;

    // Loop through each muscle in the list of muscle_names
    for (int i=0; i<numConstraints; ++i)
    {
		Constraint &thisConstraint = cs[_constraintIndex[i]];
        
		SimTK::ConstraintIndex toeConstraintIndex(0);
		SimTK::Constraint& simConstraint = _model->updMatterSubsystem().updConstraint(toeConstraintIndex);

		// number of bodies being directly constrained
		int ncb = simConstraint.getNumConstrainedBodies();
		// number of mobilities being directly constrained
		int ncm = simConstraint.getNumConstrainedU(s);

		SimTK::Vector_<SimTK::SpatialVec> bodyForcesInAncestor(ncb);
		bodyForcesInAncestor.setToZero();
		SimTK::Vector mobilityForces(ncm, 0.0);

		Array<double> values(0.0,6*ncb+ncm);

		if (!thisConstraint.isDisabled(s)) {
		
			thisConstraint.calcConstraintForces(s, bodyForcesInAncestor, mobilityForces);

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


		}

		if (abs(values[0]) > get_friction_coefficient()*values[1] && values[1] > 0) {
			TotalPenalty += (abs(values[0]) - get_friction_coefficient()*values[1]);
		}


    }

    return TotalPenalty;
}


//_____________________________________________________________________________
/** 
 * Returns the number of probe inputs in the vector returned by computeProbeInputs().
 */
int ConstraintSlipPenaltyProbe::getNumProbeInputs() const
{
    return 1;
}


//_____________________________________________________________________________
/** 
 * Provide labels for the probe values being reported.
 */
Array<string> ConstraintSlipPenaltyProbe::getProbeOutputLabels() const 
{
    Array<string> labels;
    labels.append(getName()+"_SlipPenalty");
    return labels;
}

