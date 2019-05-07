cmake_minimum_required(VERSION 3.13.0)
project(os LANGUAGES CXX CUDA)
set(OS_VERSION 0.0.0)
set(LOCAL_INCLUDE ${CMAKE_SOURCE_DIR}/include)
find_package (Threads)
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()

enable_testing()

set(lgtestpath ${CMAKE_SOURCE_DIR}/vendor/gtest/build/lib/libgtest.a)
set(cryptlib_pre ${CMAKE_SOURCE_DIR}/include/cryptly)
set(cryptlib_include ${CMAKE_SOURCE_DIR}/include)

add_library(cryptlib STATIC
	${cryptlib_pre}/cuda_ops.cu
	${cryptlib_pre}/kernels/add.cu
	${cryptlib_pre}/kernels/dot_product.cu
)

target_include_directories(cryptlib PUBLIC ${cryptlib_include})
target_compile_features(cryptlib PUBLIC cxx_std_14)
set_target_properties(cryptlib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

file(GLOB_RECURSE CUDA_TESTS ${CMAKE_SOURCE_DIR}/include/**/*spec.cxx)
message(STATUS "Found tests: ${CUDA_TESTS}")

foreach(testpath ${CUDA_TESTS})

	get_filename_component(testname ${testpath} NAME_WE)
	set(target ${testname}.spec)
	add_executable(${target} ${testpath})
	set_property(TARGET ${target}
			PROPERTY CUDA_SEPARABLE_COMPILATION ON)
	target_link_libraries(${target} PRIVATE cryptlib ${lgtestpath} ${CMAKE_THREAD_LIBS_INIT})
	target_include_directories(${target} SYSTEM
			PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/vendor/gtest/googletest/include>)
	add_test(${target} ${target})
endforeach()
message(STATUS "My ${ARCH}")
