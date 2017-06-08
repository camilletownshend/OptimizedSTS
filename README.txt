# OptimizedSTS
CMA-ES predictive optimization of sit to stand using the Ashby OpenSim model, and working off of Carmichael Ong's Standing Long Jump project program.

------------
Introduction
------------
This zip contains files needed to build an executable to perform a dynamic optimization
of a sit to stand. Details of the framework can be found in the original research which this
project is based on:

Ong CF, Hicks JL, Delp SL: Simulation-Based Design for Wearable Robotic Systems: An 
Optimization Framework for Enhancing a Standing Long Jump. IEEE Transactions on 
Biomedical Engineering, accepted, 2015.

------------
Dependencies
------------
OpenSim and Simbody must be built first. This has been tested to work with OpenSim 3.2, Simbody 3.3.1 
and OpenSim 3.3. Note that answers will change depending on your personal development environment.

----------------------------------
Provided files and other resources
----------------------------------
There are three folders given in the source package:

1. "Common" contains .h and .cpp files that have Probes, EventHandlers, and the 
   MuscleLikeCoordinateActuator that are not in the OpenSim distribution. It also 
   has files for the Covariance Matrix Adaptation (CMA) optimizer.
2. "InitialFiles" has a provided model file 
   and a controls file for a sit to stand (InitialControllerParameters.sto). 
   Use this controls file as a template to make new controls files.
3. "nominal" contains a .cpp file that creates the main() program to load a model 
   and either perform an optimization or run a forward simulation with various analyses attached.
   This file assumes that only the physiological actuators (i.e., MuscleLikeCoordinateActuators)
   are used.

The model provided is based on a model previously developed by Ashby and Delp. 
Ashby used the model to study the role of arms in a long jump in his 2006 paper, 
and a more detailed description is provided in his thesis. Differences in this 
implementation include:
1. Ligament forces are modeled using CoordinateLimitForce in OpenSim.
2. Two point constraints are used at the toe and the heel, rather than just the toe.

Since this problem is a non-convex, non-linear optimization, the CMA optimizer is
provided and used as default. For more information about how this optimizer works, 
please refer to its Wikipedia page or tutorial slides. The CMA wrapper provided 
here also has an option to enable MPI to parallelize calculations for the optimization.
1. Wikipedia: https://en.wikipedia.org/wiki/CMA-ES
2. Tutorial slides: https://www.lri.fr/~hansen/gecco2013-CMA-ES-tutorial.pdf

---------------------------
CMake Configuration Options
---------------------------
The root directory of the package has a CMakeLists.txt file that can be used to build 
the project. Key options include:

1. ENABLE_MPI: If checked, the project will be able to use the Message Passing Interface (MPI) 
   standard to allow parallelization of the CMA optimizer. 
2. MPI_INSTALL_DIR: Used for Windows OS. Assumes that the Microsoft HPC Pack has been installed 
   and will be used as the MPI implementation. (If using UNIX, this will be ignored and the MPICH 
   libraries are assumed to be used).
3. NameSpace: Leave blank if building from OpenSim source. Use "OpenSim_" when linking to pre-built OpenSim libraries.
4. OPENSIM_INSTALL_DIR: Point to the install directory for OpenSim.

--------------------
Building the project
--------------------
Building the ALL_BUILD project will build two projects:

1. Executables - StandingLongJump_nodal: This will build StandingLongJump_nodal.exe. This executable 
   is used for running both a single forward simulation or for performing an optimization.
2. Libraries - sljCommon: This will build a dynamic linked library (sljCommon.dll or sljCommon.so) 
   that contains extra classes not included in the distribution of OpenSim. This includes Probes, 
   EventHandlers, and the MuscleLikeCoordinateActuator needed for the model and optimization as well 
   as an implementation of the CMA optimizer.

-------------------------------------------------
Running a single forward simulation with analyses
-------------------------------------------------
To run a single forward simulation with analyses, we need to have:
1. The executable, StandingLongJump_nodal.exe
2. The common library, sljCommon.dll. This needs to either be in the same directory as the executable or on the system path for Windows.
3. Input model file without controls (i.e. AshbyModel_twoConstraints.osim)
4. Controls file (i.e. InitialControllerParameters.sto)

Options:
-m <modelfile.osim>: Specify a model file to use. Default is AshbyModelv8_twoConstraints.osim if not specified. 
-cf <controlsfile.sto>: Specify a controls file to use. Default is InitialControllerParameters.sto if not specified.

The executable can be run from the command prompt. As an example, you can run the executable as:
StandingLongJump_nodal.exe -m newmodel.osim -cf optimizedcontrollerparameters.sto

This will output the following files:
1. forceReporter.sto: Output from a ForceReporter.
2. fwdLog.txt: Contains information about length of simulation, integrator tolerance, 
   and values that are used to calculate the objective function.
3. PointKinematics_HEEL_*.sto: PointKinematics analysis of the heel.
4. PointKinematics_TOE_*.sto: PointKinematics analysis of the toe.

To visualize this in the OpenSim GUI:
1. Open the model "AshbyModel_twoConstraints.osim" in OpenSim.
2. Tools > Analyze
4. Select the States option and select the "states.sto" file created from the forward run of the 
executable program. It should be in the directory where you ran the executable (build/release/states.sto).
5. Check the box "Solve for equilibrium for actuator states".
6. The time range should have populated automatically (0 to 1.6).
7. Go to the Analyses tab and Add > a kinematics analysis.
8. Make sure you know where your output files are going. I made a "Runs" directory in which I made 
a folder with today's date and saved the output in that.
9. Click run. Close the Analyze tool when it's done.
10. Now you can click play and watch the simulation run!

--------------------------
Performing an optimization
--------------------------
To perform an optimization, we need to have:
1. The executable: SitToStand_nodal.exe
2. The common library: stsCommon.dll (or stsCommon.so). This needs to either be in the 
   same directory as the executable or on the system path for Windows.
3. Input model file without controls (i.e., AshbyModel_seatConstraint.osim)
4. Controls file with the initial guess for the optimizer (i.e., InitialControllerParameters.sto)

Options:
-m <modelfile.osim>: Specify a model file to use. Default is AshbyModel_twoConstraints.osim if not specified. 
-cf <controlsfile.sto>: Specify a controls file to use. Default is InitialControllerParameters.sto if not specified.
-opt: Flag that specifies that an optimization will be performed. This flag must be used when performing an optimization.
-lambda <int>: CMA optimization parameter specifying number of samples per generation. Default is 100 if not specified.
-sigma <positive real number>: CMA optimization parameter specifying constant scaling term of initial covariance matrix.
-maxIters <int>: Optimization parameter limiting number of outer iterations (a.k.a. generations for CMA). Default is 1000 if not specified.
-resume: If resumecmaes.dat file is present and the -opt flag is used, this flag specifies the CMA optimizer to resume from a previous 
         optimization corresponding to the file. This will override any -sigma setting, but -lambda and -maxIters can still be specified. 
         If resumecmaes.dat is not present, then the optimizer will begin a new optimization. If the -opt flag is not used, then the -resume 
         flag will have no effect and a single forward simulation will be performed.

The executable can be run from the command prompt. As an example, you can run the executable as:
"StandingLongJump_nodal.exe -opt -m newmodel.osim -cf initialcontrollerparameters.sto -lambda 50 -sigma 0.01 -maxIters 2000"

If the ENABLE_MPI flag was set to true in CMake and the proper MPI libraries are linked, you can use more than one thread to 
parallelize calculations. An example way to run this (in this example using 8 threads) is:
"mpiexec -n 8 SitToStand_nodal.exe -opt -m newmodel.osim -cf initialcontrollerparameters.sto -lambda 50 -sigma 0.01 -maxIters 2000"

An example to resume a previous CMA optimization is to run the following in the command prompt:
"mpiexec -n 8 SitToStand_nodal.exe -opt -resume -m newmodel.osim -cf controllerparameters.sto -lambda 50 -maxIters 2000"

This will output the following files specific to this project:
1. ControllerParameters_gen<int>.sto: This is the controller parameters corresponding to the best objective function value for the 
   generation number. This is set to print the 1st and every 100th generations.
2. optimizedControls.sto: Controller parameters corresponding to the best objective function found.
3. optimizedModel.osim: Model with the controller and control parameters corresponding to optimizedControls.sto.
4. optLog.txt: Lists the options and parameters used for the optimization.
5. outputlog.txt: Tab-delimited file with information about the progression of the optimizer. The columns contain: 
   1) generation number, 2) best objective function value, 3) sigma for the current generation, and 
   4) corresponding controls for the best objective function. This is printed for the 1st and every 50th generations.
6. resumecmaes.dat: Generated by the CMA optimizer to allow the user to resume an optimization.