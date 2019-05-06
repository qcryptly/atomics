cmake_minimum_required(VERSION 3.13.0)
project(os LANGUAGES CXX CUDA C)
set(OS_VERSION 0.0.0)
set(LOCAL_INCLUDE ${CMAKE_SOURCE_DIR}/include)
find_package (Threads)
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()

message(STATUS "LOCAL_INCLUDE=${LOCAL_INCLUDE}")
enable_testing()

#file(GLOB SRCS ${CMAKE_SOURCE_DIR}/src/**/*.cxx)
#message(STATUS "Found srcs: ${SRCS}")
#foreach(srcpath ${SRCS})
#  get_filename_component(srcname ${srcpath} NAME_WE)
#  set(target run-${srcname})
#  add_executable(${target} ${srcpath})
#  target_compile_options(${target} PRIVATE -g -Wall -std=c++2a)
#  target_include_directories(${target} SYSTEM
#                             PUBLIC $<BUILD_INTERFACE:${LOCAL_INCLUDE}>)
#endforeach()
#
set(lgtestpath ${CMAKE_SOURCE_DIR}/vendor/gtest/build/lib/libgtest.a)
set(cryptlib_pre ${CMAKE_SOURCE_DIR}/include/cryptly)
set(cryptlib_include ${CMAKE_SOURCE_DIR}/include)

add_library(cryptlib SHARED
	${cryptlib_pre}/cuda_ops.cu
	${cryptlib_pre}/kernels/add.cu
	${cryptlib_pre}/kernels/dot_product.cu
)
target_include_directories(cryptlib PUBLIC ${cryptlib_include})
target_compile_features(cryptlib PUBLIC cxx_std_14)
set_target_properties(cryptlib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(cuda_ops_test ${cryptlib_pre}/cuda_ops.spec.cxx)
set_property(TARGET cuda_ops_test
		PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(cuda_ops_test PRIVATE cryptlib ${lgtestpath} ${CMAKE_THREAD_LIBS_INIT})
target_include_directories( cuda_ops_test SYSTEM
		PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/vendor/gtest/googletest/include>)

#file(GLOB CUDA_TESTS ${CMAKE_SOURCE_DIR}/include/**/*spec.cxx)
#message(STATUS "Found tests: ${TESTS}")
#
#foreach(testpath ${TESTS})
#  get_filename_component(testname ${testpath} NAME_WE)
#  set(target ${testname}.spec)
#  add_executable(${target} ${testpath})
#  target_compile_options(${target} PRIVATE -g -Wall -pthread -std=c++2a)
#  target_link_libraries(${target} PUBLIC ${lgtestpath} ${CMAKE_THREAD_LIBS_INIT})
#  target_include_directories(
#    ${target} PUBLIC
#    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
#  target_include_directories(
#    ${target} SYSTEM
#    PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/vendor/gtest/googletest/include>)
#  add_test(${target} ${target})
#  message(STATUS "Building unit test ${target}")
#endforeach()
#message(STATUS "My ${ARCH}")
