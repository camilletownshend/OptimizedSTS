#ifndef OPENSIM_VERTICAL_GRF_PROBE_H_
#define OPENSIM_VERTICAL_GRF_PROBE_H_
/* -------------------------------------------------------------------------- *
 *                    OpenSim:  VerticalGRFProbe.h                            *
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

//============================================================================
// INCLUDES
//============================================================================
#include "sljCommonDLL.h"
#include <OpenSim/Simulation/Model/Probe.h>


namespace OpenSim {

//=============================================================================
//=============================================================================
/**
 * VerticalGRFProbe is a ModelComponent Probe for computing the 
 * Vertical GRF from PointConstraints
 *
 * @author Carmichael Ong
 */

class Model;

class SLJCOMMON_API VerticalGRFProbe : public Probe
{
    OpenSim_DECLARE_CONCRETE_OBJECT(VerticalGRFProbe, Probe);

public:
//==============================================================================
// PROPERTIES
//==============================================================================
    /** @name Property declarations
    These are the serializable properties associated with this class. **/
    /**@{**/
    /** List of VerticalGRFs to probe.  **/
    OpenSim_DECLARE_LIST_PROPERTY(constraint_names, std::string,
        "Specify a list of model VerticalGRFs whose activation should be calculated."
        "Use 'all' to probe all coordinate limit forces.");

    /** Element-wise power exponent to apply to each coordinate limit force prior to the Probe operation. 
    For example, if two coordinate limit forces M1 and M2 are given in coordinate_liimit_force_names, then the
    Probe value will be equal to Force_M1^exponent + Force_M2^exponent.  **/
    OpenSim_DECLARE_PROPERTY(exponent, double,
        "Element-wise power exponent to apply to each coordinate limit force prior to the Probe operation.");
    /**@}**/

//==============================================================================
// PUBLIC METHODS
//==============================================================================
    //--------------------------------------------------------------------------
    // Constructor(s) and Setup
    //--------------------------------------------------------------------------
    /** Default constructor */
    VerticalGRFProbe();

    // Uses default (compiler-generated) destructor, copy constructor, and copy
    // assignment operator.


    //--------------------------------------------------------------------------
    // Computation
    //--------------------------------------------------------------------------
    /** Compute the Coordinate Limit Force */
    SimTK::Vector computeProbeInputs(const SimTK::State& state) const OVERRIDE_11;

    /** Returns the number of probe inputs in the vector returned by computeProbeInputs(). */
    int getNumProbeInputs() const OVERRIDE_11;

    /** Returns the column labels of the probe values for reporting. 
        Currently uses the Probe name as the column label, so be sure
        to name your probe appropiately! */
    Array<std::string> getProbeOutputLabels() const OVERRIDE_11;


//==============================================================================
// PRIVATE
//==============================================================================
private:
    // The index inside OpenSim::ConstraintSet that corresponds to each Actuator
    // being probed.
    SimTK::Array_<int> _constraintIndex;

    //--------------------------------------------------------------------------
    // ModelComponent Interface
    //--------------------------------------------------------------------------
    void connectToModel(Model& aModel) OVERRIDE_11;
    
    void setNull();
    void constructProperties();

//==============================================================================
};	// END of class VerticalGRFProbe
//==============================================================================
//==============================================================================

} // end of namespace OpenSim

#endif // OPENSIM_COORDINATE_LIMIT_FORCE_PROBE_H_


