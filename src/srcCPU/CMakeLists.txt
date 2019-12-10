
include_directories(${THRUST_INCLUDE_DIR}) #specifies path that compiler looks in while searching for header files
 
include_directories(${COMMON_INCLUDES})
include_directories(${MATLAB_DIR}/extern/include)
include_directories(${MATLAB_DIR}/bin/glnxa64)
include_directories(${MATLAB_DIR}/extern/include/MatlabEngine)
add_library(SceInit STATIC
        CellInitHelper.h
        CellInitHelper.cpp
        Signal2D.h
		Signal2D.cpp
		Signal_Calculator.hpp
		Signal_Calculator.cpp
		Tissue.hpp
		Tissue.cpp
		Mesh.hpp
		Mesh.cpp
		Cell.hpp
		Cell.cpp
		common.cpp
		common.hpp
        ${COMMON_INCLUDES}/Point2D.cpp
        ${COMMON_INCLUDES}/Point2D.h
        ${COMMON_INCLUDES}/ConfigParser.h
        ${COMMON_INCLUDES}/ConfigParser.cpp
        ${COMMON_INCLUDES}/GeoVector.h
        ${COMMON_INCLUDES}/GeoVector.cpp
)

target_compile_features(SceInit PUBLIC cxx_std_11)

target_link_libraries(SceInit
			${MATLAB_DIR}/bin/glnxa64/libeng.so
			${MATLAB_DIR}/bin/glnxa64/libmx.so 
			${MATLAB_DIR}/bin/glnxa64/libmex.so
			${MATLAB_DIR}/bin/glnxa64/libmat.so
			${MATLAB_DIR}/extern/bin/glnxa64/libMatlabDataArray.so
			${MATLAB_DIR}/extern/bin/glnxa64/libMatlabEngine.so)