#MESSAGE(${CGAL_INCLUDE_DIRS})
#include_directories(${CGAL_INCLUDE_DIRS})
include_directories(${THRUST_INCLUDE_DIR}) #specifies path that compiler looks in while searching for header files
 
include_directories(${COMMON_INCLUDES})
include_directories(/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018a/extern/include)
include_directories(/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018a/bin/glnxa64)
include_directories(/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018a/extern/include/MatlabEngine)
add_library(SceInit STATIC
        CellInitHelper.h
        CellInitHelper.cpp
        Signal2D.h
		Signal2D.cpp
        ${COMMON_INCLUDES}/Point2D.cpp
        ${COMMON_INCLUDES}/Point2D.h
        ${COMMON_INCLUDES}/ConfigParser.h
        ${COMMON_INCLUDES}/ConfigParser.cpp
        ${COMMON_INCLUDES}/GeoVector.h
        ${COMMON_INCLUDES}/GeoVector.cpp
)

target_compile_features(SceInit PUBLIC cxx_std_11)

target_link_libraries(SceInit
			/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018a/bin/glnxa64/libeng.so
			/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018a/bin/glnxa64/libmx.so 
			/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018a/bin/glnxa64/libmex.so
			/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018a/bin/glnxa64/libmat.so
			/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018a/extern/bin/glnxa64/libMatlabDataArray.so
			/opt/linux/centos/7.x/x86_64/pkgs/MATLAB/R2018a/extern/bin/glnxa64/libMatlabEngine.so)