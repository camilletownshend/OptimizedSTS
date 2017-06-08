#ifndef _SimTK_CMA_OPTIMIZER_H_
#define _SimTK_CMA_OPTIMIZER_H_

#include "sljCommonDLL.h"
#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/OptimizerRep.h" 
#include "cmaes_interface.h"

namespace SimTK {


class SLJCOMMON_API CMAOptimizer : public Optimizer::OptimizerRep {
public:
    ~CMAOptimizer(); 

    CMAOptimizer(const OptimizerSystem& sys); 

    Real optimize(  SimTK::Vector &results );
    OptimizerRep* clone() const;

	// override for base class in Simbody 3.4+
	OptimizerAlgorithm getAlgorithm() const {
		return OptimizerAlgorithm(5);
	}

private:
	cmaes_t evo;
	int lambda; // number of samples per iteration
	double sigma; // initial stepsize
	bool resume;  
	bool enableMPI;  

	void slave(); 
	Real master( SimTK::Vector &results ); 
};

} // namespace SimTK

#endif //_SimTK_CMA_OPTIMIZER_H_

