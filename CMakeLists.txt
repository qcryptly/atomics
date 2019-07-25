cmake_minimum_required(VERSION 3.13.0)
project(os LANGUAGES CXX CUDA)
set(OS_VERSION 0.0.0)
set(LOCAL_INCLUDE ${CMAKE_SOURCE_DIR}/include)
find_package (Threads)
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()

enable_testing()

find_package(CUDA 10.1 REQUIRED)
find_package(Boost REQUIRED COMPONENTS system python3)
find_package(PythonLibs ${PY_VERSION} REQUIRED)

set(lgtestpath ${CMAKE_SOURCE_DIR}/vendor/gtest/build/lib/libgtest.a)
set(wrapperlib_include ${CMAKE_SOURCE_DIR}/include)

file(GLOB_RECURSE WRAPPER_FILES ${wrapperlib_include}/**/*.cu)

add_library(wrapperlib SHARED ${WRAPPER_FILES})
target_include_directories(wrapperlib PUBLIC ${wrapperlib_include})
set_target_properties(wrapperlib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

file(GLOB_RECURSE CUDA_TESTS ${CMAKE_SOURCE_DIR}/include/**/*spec.cxx)
message(STATUS "Found tests: ${CUDA_TESTS}")

foreach(testpath ${CUDA_TESTS})
	get_filename_component(testname ${testpath} NAME_WE)
	set(target ${testname}.spec)
	add_executable(${target} ${testpath})
	set_property(TARGET ${target}
			PROPERTY CUDA_SEPARABLE_COMPILATION ON)
	target_include_directories(${target} SYSTEM
			PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/vendor/gtest/googletest/include>
			${Boost_INCLUDE_DIRS})
	target_link_libraries(${target} PRIVATE wrapperlib ${lgtestpath} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} ${CUDA_LIBRARIES})
	add_test(${target} ${target})
endforeach()
message(STATUS "My ${ARCH}")

file(GLOB_RECURSE PYTHON_SRCS ${CMAKE_SOURCE_DIR}/pythonsrc/*.cxx)
message(STATUS "Found python sources")

# Create python cryptmat library
add_library(cryptmat SHARED ${PYTHON_SRCS})
set_target_properties(cryptmat PROPERTIES PREFIX "")
set_property(TARGET cryptmat PROPERTY CUDA_SEPERABLE_COMPILATION ON)
target_include_directories(cryptmat SYSTEM 
  PUBLIC ${Boost_INCLUDE_DIRS} /usr/include/python3.7m)
target_link_libraries(cryptmat PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} boost_numpy3 wrapperlib)
configure_file(${CMAKE_SOURCE_DIR}/pythonsrc/__init__.py ${CMAKE_BINARY_DIR}/__init__.py COPYONLY)
