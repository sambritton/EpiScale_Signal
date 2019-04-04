include_directories(${GOOGLETEST_DIR}/include/ ${COMMON_INCLUDES} ${GPU_SOURCE_DIR} ${CPU_SOURCE_DIR})
include_directories(${THRUST_INCLUDE_DIR}) #specifies path that compiler looks in while searching for header files

include_directories(${Matlab_INCLUDE_DIRS})

include_directories(${COMMON_CUDA_INCLUDES}) #specifies path that compiler looks in while searching for header files

set(DISK_PROJECT_NAME runDiscSimulation)
MESSAGE("Now building " ${DISK_PROJECT_NAME})
file(GLOB SRC_FILES ./discMain.cpp)
MESSAGE("Compiling these files: " ${SRC_FILES})

add_executable(${DISK_PROJECT_NAME} ${SRC_FILES})
add_dependencies(${DISK_PROJECT_NAME} gtest)
add_dependencies(${DISK_PROJECT_NAME} gtest_main)
add_dependencies(${DISK_PROJECT_NAME} SceGPU)

if(NOT WIN32)
        #target_link_libraries(${DISK_PROJECT_NAME} pthread)
        target_link_libraries(${DISK_PROJECT_NAME} SceGPU)
        target_link_libraries(${DISK_PROJECT_NAME} SceInit)
        target_link_libraries(${DISK_PROJECT_NAME} SceUtils)
        #target_link_libraries(${DISK_PROJECT_NAME} SceUtilsGPU)

        target_link_libraries(${DISK_PROJECT_NAME}
                ${GOOGLETEST_DIR}/libgtest.a
                ${GOOGLETEST_DIR}/libgtest_main.a
        )

	#	target_link_libraries(${DISK_PROJECT_NAME}
	#				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libeng.so
	#				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmx.so 
	#				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmex.so
	#				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmat.so )
endif()

set(MODEL_TEST_NAME runModelTest)
MESSAGE("Now building " ${MODEL_TEST_NAME})
file(GLOB SRC_FILES ./modelTest.cpp)
#MESSAGE("Compiling these files: "${SRC_FILES})
add_executable(${MODEL_TEST_NAME} ${SRC_FILES})
add_dependencies(${MODEL_TEST_NAME} gtest)
add_dependencies(${MODEL_TEST_NAME} gtest_main)
add_dependencies(${MODEL_TEST_NAME} SceGPU)

if(NOT WIN32)
        #target_link_libraries(${MODEL_TEST_NAME} pthread)
        target_link_libraries(${MODEL_TEST_NAME} SceGPU)
        target_link_libraries(${MODEL_TEST_NAME} SceInit)
        target_link_libraries(${MODEL_TEST_NAME} SceUtils)
        #target_link_libraries(${MODEL_TEST_NAME} SceUtilsGPU)

        target_link_libraries(${MODEL_TEST_NAME}
                ${GOOGLETEST_DIR}/libgtest.a
                 ${GOOGLETEST_DIR}/libgtest_main.a)

		
endif()

set(LASER_PROJECT_NAME LaserAblationTest)
MESSAGE("Now building " ${LASER_PROJECT_NAME})
file(GLOB SRC_FILES ./laserAblationTest.cpp)
#MESSAGE("Compiling these files: "${SRC_FILES})
add_executable(${LASER_PROJECT_NAME} ${SRC_FILES})
add_dependencies(${LASER_PROJECT_NAME} gtest)
add_dependencies(${LASER_PROJECT_NAME} gtest_main)
add_dependencies(${LASER_PROJECT_NAME} SceGPU)

if(NOT WIN32)
        target_link_libraries(${LASER_PROJECT_NAME} SceGPU)
        target_link_libraries(${LASER_PROJECT_NAME} SceInit)
        target_link_libraries(${LASER_PROJECT_NAME} SceUtils)

        target_link_libraries(${LASER_PROJECT_NAME}
                ${GOOGLETEST_DIR}/libgtest.a
                ${GOOGLETEST_DIR}/libgtest_main.a
        )

endif()

set(DISK_M_PROJECT_NAME runDiscSimulation_M)
MESSAGE("Now building " ${DISK_M_PROJECT_NAME})
file(GLOB SRC_FILES ./discMain_M.cu)
#MESSAGE("Compiling these files: "${SRC_FILES})
add_executable(${DISK_M_PROJECT_NAME} ${SRC_FILES})
add_dependencies(${DISK_M_PROJECT_NAME} gtest)
add_dependencies(${DISK_M_PROJECT_NAME} gtest_main)
add_dependencies(${DISK_M_PROJECT_NAME} SceGPU)

if(NOT WIN32)

        #target_link_libraries(${DISK_M_PROJECT_NAME} pthread)
        target_link_libraries(${DISK_M_PROJECT_NAME} SceGPU)
        target_link_libraries(${DISK_M_PROJECT_NAME} SceInit)
        target_link_libraries(${DISK_M_PROJECT_NAME} SceUtils)
        #target_link_libraries(${DISK_M_PROJECT_NAME} SceUtilsGPU)

        target_link_libraries(${DISK_M_PROJECT_NAME}
                ${GOOGLETEST_DIR}/libgtest.a
                ${GOOGLETEST_DIR}/libgtest_main.a
        )

	target_link_libraries(${DISK_M_PROJECT_NAME}
				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libeng.so
				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmx.so 
				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmex.so
				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmat.so
				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/bin/glnxa64/libMatlabDataArray.so
				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/bin/glnxa64/libMatlabEngine.so)
							
	#	target_link_libraries(${DISK_M_PROJECT_NAME}
	#				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libeng.so
	#				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmx.so 
	#				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmex.so
	#				/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/bin/glnxa64/libmat.so )
endif()

