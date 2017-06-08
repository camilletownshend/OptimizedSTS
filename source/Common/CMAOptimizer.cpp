#include "CMAOptimizer.h"
#include <iostream>
#include <fstream>
#ifdef USE_MPI
#include <mpi.h>
#endif
#define DIETTAG 2

using std::cout;
using std::endl;

namespace SimTK {

Optimizer::OptimizerRep* CMAOptimizer::clone() const {
    return( new CMAOptimizer(*this) );
}


CMAOptimizer::CMAOptimizer( const OptimizerSystem& sys )
    : OptimizerRep( sys ), lambda(0), sigma(0.3), resume(false), enableMPI(false) 
{
     /* internal flags for CMA */

     if( sys.getNumParameters() < 1 ) {
        const char* where = "Optimizer Initialization";
        const char* szName = "dimension";
        SimTK_THROW5(SimTK::Exception::ValueOutOfRange, szName, 1,  sys.getNumParameters(), INT_MAX, where); 
     }

}
CMAOptimizer::~CMAOptimizer() {
}

void CMAOptimizer::slave() {
#ifdef USE_MPI
	MPI_Status status;
    const OptimizerSystem& sys = getOptimizerSystem();
	int nvars = sys.getNumParameters();
	Vector controls; 
	controls.resize(nvars); 
	double* msg = new double[nvars+1];

	while (1) {
		MPI_Recv(msg,
				nvars+1 ,
				MPI_DOUBLE, 0,
				MPI_ANY_TAG,
				MPI_COMM_WORLD, &status);

		if (status.MPI_TAG == DIETTAG) {
			delete msg;
			MPI_Finalize(); 
			exit(0); 
		}

		for (int i = 0; i < nvars; i++) {
			controls[i] = msg[i];
		}
		double res; 
		sys.objectiveFunc(controls, true, res);

		double result[2];
		result[0] = res;
		result[1] = msg[nvars];

		MPI_Send(result, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
#endif
}

Real CMAOptimizer::master( Vector &results ) {
	std::ofstream outputFile("outputlog.txt");
    Real f = HUGE_VAL; 
    const OptimizerSystem& sys = getOptimizerSystem();
    int n = sys.getNumParameters();

    Real *lowerLimits, *upperLimits;

    if( sys.getHasLimits() ) {
        sys.getParameterLimits( &lowerLimits, &upperLimits );
    }

	int numsamples = 0;
    getAdvancedIntOption("lambda", numsamples );
	if (numsamples == 0) {
		numsamples = 4+int(3*std::log(double(n))); 
	}
	
	double stepsize = 0;
    getAdvancedRealOption("sigma", stepsize ); 
	if (stepsize == 0.0) {
		stepsize = 0.1; 
	}
	double* stepsizeArray = new double[n];
	for (int i = 0; i < n; i++) {
		stepsizeArray[i] = stepsize;  
	}
	
	bool isresume = false; 
    getAdvancedBoolOption("resume", isresume );
	
	bool usempi = false; 
    getAdvancedBoolOption("enableMPI", usempi );
	 
	double* funvals = cmaes_init(&evo, n, &results[0], stepsizeArray, 0, numsamples, "non"); 
	if (isresume) {
		cmaes_resume_distribution(&evo, "resumecmaes.dat"); 
	}

	evo.sp.stopMaxIter = maxIterations;
	std::vector<double*> msgs; 
	msgs.resize(numsamples);
	for (unsigned int i = 0; i < msgs.size(); i++) {
		msgs[i] = new double[n+1];
	} 
	
std::cout << "maxIterations " << evo.sp.stopMaxIter << std::endl;
std::cout << "step size " << cmaes_Get(&evo, "sigma") << std::endl;
std::cout << "lambda " << numsamples << std::endl;
	
	while (!cmaes_TestForTermination(&evo)) {
		double*const* pop = cmaes_SamplePopulation(&evo);
    	if( sys.getHasLimits() ) {
			for (int i = 0; i < cmaes_Get(&evo, "popsize"); i++) {
				bool feasible = false; 
				while (!feasible) {
					feasible = true; 
					for (int j = 0; j < n; j++) {
						if (pop[i][j] < lowerLimits[j] || 
								pop[i][j] > upperLimits[j]) {
							feasible = false; 
							pop = cmaes_ReSampleSingle(&evo, i); 
							break; 
						}
					}
				}
			}
		}

#ifdef USE_MPI
		if (usempi) {
			int nprocs; 	
			MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

			double result[2];
			MPI_Status status;
			int numToReceive = msgs.size();
			unsigned int nextToSend = 0; 
			for (int i = 0; i < msgs.size(); i++) {
				for (int j = 0; j < n; j++) {
					msgs[i][j] = pop[i][j]; 
				}
				msgs[i][n] = i; 
			}


			for (int i = 1; i < nprocs; i++) {
				MPI_Send(msgs[nextToSend++],
						n+1,
						MPI_DOUBLE,
						i,
						1,
						MPI_COMM_WORLD);
				if (nextToSend >= msgs.size())
					break;
			}

			while (nextToSend < msgs.size()) {
				MPI_Recv(result,
						2,
						MPI_DOUBLE,
						MPI_ANY_SOURCE,
						MPI_ANY_TAG,
						MPI_COMM_WORLD,
						&status);

				MPI_Send(msgs[nextToSend++],
						n+1,
						MPI_DOUBLE,
						status.MPI_SOURCE,
						1,
						MPI_COMM_WORLD);

				funvals[int(result[1])] = result[0];
				numToReceive--;
			}

			while (numToReceive > 0) {
				MPI_Recv(result,
						2,
						MPI_DOUBLE,
						MPI_ANY_SOURCE,
						MPI_ANY_TAG,
						MPI_COMM_WORLD,
						&status);
				funvals[int(result[1])] = result[0];
				numToReceive--;
			}
		} else {
			for (int i = 0; i < cmaes_Get(&evo, "lambda"); i++) {
				objectiveFuncWrapper(n, pop[i], true, &funvals[i], this);
			}
		}
#else
			for (int i = 0; i < cmaes_Get(&evo, "lambda"); i++) {
				objectiveFuncWrapper(n, pop[i], true, &funvals[i], this);
			}
#endif
		
		cmaes_UpdateDistribution(&evo, funvals);
		
		const double* optx = cmaes_GetPtr(&evo, "xbestever");
		for (int i = 0; i < n; i++) {
			results[i] = optx[i]; 
		}
		f = cmaes_Get(&evo, "fbestever");

        /* Printing out internal information from optimizer to check on progress */
		if (int(evo.gen) % 50 == 0 || int(evo.gen) == 1) {
			std::cout << "best val: " << f << " (gen count = " << evo.gen << ") (sigma = " << cmaes_Get(&evo, "sigma") << ")" << std::endl;
			outputFile << evo.gen << "\t" << f << "\t" << cmaes_Get(&evo, "sigma") << "\t" << results << std::endl;
		}

		if (int(evo.gen) % 100 == 0 || int(evo.gen) == 1) {
			cmaes_WriteToFile(&evo, "resume", "resumecmaes.dat");

			int numActuators = 4;
			std::string fileName = "ControllerParameters_gen";
			fileName += String(static_cast<int>(evo.gen));
			fileName += ".sto";
			std::ofstream currentOutputFile(fileName.c_str(), std::ofstream::out);

			// HEADERS
			currentOutputFile << "ControllerParameters_gen" << evo.gen << endl;
			currentOutputFile << "version=1" << endl;
			currentOutputFile << "nRows=" << results.size()/numActuators << endl;
			currentOutputFile << "nColumns=" << numActuators+1 << endl;
			currentOutputFile << "inDegrees=no" << endl;
			currentOutputFile << "endheader" << endl;

			currentOutputFile << "time" << "\t";
			currentOutputFile << "ankle_MLCA" << "\t";
			currentOutputFile << "knee_MLCA" << "\t";
			currentOutputFile << "hip_MLCA" << "\t";
			currentOutputFile << "shoulder_MLCA" << "\t";
			currentOutputFile << endl;
			currentOutputFile.precision(21);

			double timeStep = 0.05;
			
			for (int i = 0; i < 132 / numActuators; i++) {
				currentOutputFile << timeStep*i << "\t";
				for (int j = 0; j < numActuators; j++) {
					currentOutputFile << results[i + 132 / numActuators] << "\t";
				}
				currentOutputFile << endl;
			}
			currentOutputFile.close();

		}
	}
	cmaes_WriteToFile(&evo, "resume", "resumecmaes.dat");

	for (unsigned int i = 0; i < msgs.size(); i++) {
		delete msgs[i];
	}
	delete stepsizeArray;
	
	return f;  
}

Real CMAOptimizer::optimize( Vector &results ) {
#ifdef USE_MPI
	int myrank = 0; 
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	if (myrank > 0) {
		slave(); 
		return 0; 
	}
	else {
		return master(results); 
	}
#else
	return master(results); 
#endif
}

} // namespace SimTK
