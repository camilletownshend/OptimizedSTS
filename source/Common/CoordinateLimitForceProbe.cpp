/* -------------------------------------------------------------------------- *
 *                   OpenSim:  CoordinateLimitForceProbe.cpp                  *
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
#include "CoordinateLimitForceProbe.h"
#include <OpenSim/Simulation/Model/Force.h>
#include <OpenSim/Simulation/Model/ForceSet.h>
#include <algorithm>

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
CoordinateLimitForceProbe::CoordinateLimitForceProbe() 
{
    setNull();
    constructProperties();
}


//_____________________________________________________________________________
/**
 * Set the data members of this TaskTermProbe to their null values.
 */
void CoordinateLimitForceProbe::setNull(void)
{
    setAuthors("Carmichael Ong");
    _forceIndex.clear();
}

//_____________________________________________________________________________
/**
 * Connect properties to local pointers.
 */
void CoordinateLimitForceProbe::constructProperties(void)
{
    constructProperty_coordinate_limit_force_names();
    constructProperty_exponent(2.0);
	constructProperty_force_threshold(0.0);
}



//=============================================================================
// MODEL COMPONENT METHODS
//=============================================================================
//_____________________________________________________________________________
/**
 * Perform some set up functions that happen after the
 * object has been deserialized or copied.
 *
 * @param model OpenSim model containing this CoordinateLimitForceProbe.
 */
void CoordinateLimitForceProbe::connectToModel(Model& model)
{
    Super::connectToModel(model);
    _forceIndex.clear();

	if (getProperty_coordinate_limit_force_names().size() == 0)
		return;

    // Check to see if 'all' coordinate limit forces are selected for probing.
	if (get_coordinate_limit_force_names(0) == "all") {
        Array<string> allCoordinateLimitForceNames, allForceNames;
		const ForceSet forceSet = _model->getForceSet();
		forceSet.getNames(allForceNames);
		
		int numForces = allForceNames.getSize();
		for (int i=0; i<numForces; i++) {
			//if ( forceSet.get(allForceNames[i]).isA("CoordinateLimitForce") ) {
			if ( strcmp(forceSet.get(allForceNames[i]).getConcreteClassName().c_str(),"CoordinateLimitForce") == 0 ) {
				allCoordinateLimitForceNames.append(forceSet.get(allForceNames[i]).getName());
			} else {
				//cout << "i = " << i << ": " << forceSet.get(allForceNames[i]).getConcreteClassName() << endl;;
			}
		}
		
        set_coordinate_limit_force_names(allCoordinateLimitForceNames);
        //cout << "Set to all coordinate limit forces: " << allCoordinateLimitForceNames << endl;
    }

    // Check that each CoordinateLimitForce in the coordinate_limit_force array exists in the model.
    // If so, populate the _forceIndex internal array.
    const int numCoordinateLimitForces = getProperty_coordinate_limit_force_names().size();
    for (int i=0; i<numCoordinateLimitForces; i++) {
        const string& coordinateLimitForceName = get_coordinate_limit_force_names(i);
        const int k = model.getForceSet().getIndex(coordinateLimitForceName);
        if (k<0) {
            string errorMessage = getConcreteClassName() + ": Invalid CoordinateLimitForce '" 
                + coordinateLimitForceName + "' specified in <coordinate_limit_force_names>.";
            throw (Exception(errorMessage.c_str()));
        }
        else
            _forceIndex.push_back(k);
    }

    // Sanity check.
    //if (numCoordinateLimitForces != int(_actuatorIndex.size()))
    //    throw (Exception("Size of _actuatorIndex does not match number of actuators."));
}





//=============================================================================
// COMPUTATION
//=============================================================================
//_____________________________________________________________________________
/**
 * Compute the summed coordinate limit forces.
 */
Vector CoordinateLimitForceProbe::computeProbeInputs(const State& s) const
{
    const int numCoordinateLimitForces = getProperty_coordinate_limit_force_names().size();
    const ForceSet& fs = _model->getForceSet();
    SimTK::Vector TotalF(1);
    TotalF = 0;

    // Loop through each coordinate limit force in coordinate_limit_force_names
    for (int i=0; i<numCoordinateLimitForces; ++i)
    {
        // Get the muscle activation and append to output.
		Array<double> currentForceValues = fs[_forceIndex[i]].getRecordValues(s);
		const double force = std::max(abs(currentForceValues[0]) - get_force_threshold(), 0.0);
        TotalF(0) += std::pow(force, get_exponent());
    }

    return TotalF;
}


//_____________________________________________________________________________
/** 
 * Returns the number of probe inputs in the vector returned by computeProbeInputs().
 */
int CoordinateLimitForceProbe::getNumProbeInputs() const
{
    return 1;
}


//_____________________________________________________________________________
/** 
 * Provide labels for the probe values being reported.
 */
Array<string> CoordinateLimitForceProbe::getProbeOutputLabels() const 
{
    Array<string> labels;
    labels.append(getName()+"_Summed");
    return labels;
}

