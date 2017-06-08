SET ROOT=C:\Users\Camille\OpenSim\ME485ProjectWorkSpace\STS
SET SRC=%ROOT%\source
SET BUILD=%ROOT%\build\Release
%BUILD%\predictive_STS_nominal.exe -m %ROOT%\STS_SeatConstraint.osim -cf %ROOT%\InitialControllerParameters.sto -opt -lambda 5 -sigma 0.1 -maxIters 200