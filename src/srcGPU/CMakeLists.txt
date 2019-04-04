#MESSAGE(${GPU_SOURCE_DIR})
include_directories(${COMMON_CUDA_INCLUDES} )
MESSAGE("in scegpu COMMON_CUDA_INCLUDES: " ${COMMON_CUDA_INCLUDES})

include_directories(${COMMON_INCLUDES} ${GOOGLETEST_DIR}/include/)
include_directories(${THRUST_INCLUDE_DIR}) #specifies path that compiler looks in while searching for header files
include_directories(${CUDA_TOOLKIT_ROOT_DIR}) 

MESSAGE("in scegpu " ${THRUST_INCLUDE_DIR}/)
MESSAGE("in scegpu " ${GOOGLETEST_DIR}/include/)

add_library(SceGPU STATIC
	SceNodes.cu
	SceNodes.h
	SceCells.cu
	SceCells.h
	SimulationDomainGPU.cu
	SimulationDomainGPU.h
	Structures_Misc.h
)


target_compile_features(SceGPU PUBLIC cxx_std_11) 