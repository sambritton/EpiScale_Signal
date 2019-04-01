#MESSAGE(${CGAL_INCLUDE_DIRS})
#include_directories(${CGAL_INCLUDE_DIRS})
include_directories(${THRUST_INCLUDE_DIR}) #specifies path that compiler looks in while searching for header files
 
include_directories(${COMMON_INCLUDES})
include_directories(/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/include)
include_directories(/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64)
include_directories(/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/include/MatlabEngine)
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
			/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libeng.so
			/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmx.so 
			/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmex.so
			/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmat.so
			/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/bin/glnxa64/libMatlabDataArray.so
			/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/bin/glnxa64/libMatlabEngine.so)