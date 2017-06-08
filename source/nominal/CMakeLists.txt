SET(APPNAME predictive_STS_nominal)

# CMakeLists file for ForwardJumpingOptimization_toeConstraint\CMA

# Identify the cpp file(s) that were to be built
FILE(GLOB SOURCE_FILES *.h *.cpp *.c)
SET(SOURCE ${SOURCE_FILES})

LINK_LIBRARIES(sljCommon)

ADD_EXECUTABLE(${APPNAME} ${SOURCE})
SET_TARGET_PROPERTIES(${APPNAME} PROPERTIES PROJECT_LABEL "Executables - ${APPNAME}")