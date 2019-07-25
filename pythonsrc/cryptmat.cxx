#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <utility>
#include <iostream>
#include <tuple>
#include <unordered_map>
#include <cassert>

// Used for primitive_types we are generalizing
#include "cryptly/cuda_ops.hxx"
#include "cryptly/kernels/dot_product.hxx"

using type_lists = std::tuple<primitive_types>;

namespace py = boost::python;
namespace np = boost::python::numpy;

template<auto... index>
constexpr auto np_make_types(std::index_sequence<index...>){
    return std::make_tuple(
      (np::dtype::get_builtin<
        typename std::tuple_element<index, type_lists>::type>())...);
};

const char* greet(){return "Hello world!";}

//np::ndarray data_ex1 = np::from_data(data,dt, shape,stride,own);

template<class T>
py::object dot_impl(std::uint64_t size, void * data_a, void * data_b) {
	constexpr auto block_size = 256;
	constexpr auto grid_size = 256;
	T results;
	cryptly::kernel_ps_t kernel{block_size, grid_size};
	(cryptly::dot_product{})(kernel, size,
			reinterpret_cast<T*>(data_a),
			reinterpret_cast<T*>(data_b),
			&results);
	return py::object(results);
};

py::object dot(const np::ndarray& a, const np::ndarray& b) {
  assert(a.get_shape() == b.get_shape());
	assert(a.get_dtype() == b.get_dtype());
	assert(a.get_strides() == b.get_strides());
	
  std::uint64_t ndim = py::extract<std::uint64_t>(a.attr("ndim"));
  auto shape = a.get_shape();
	auto strides = a.get_strides();
	auto dtype_py = a.get_dtype(); 
  std::string dtype{py::extract<const char *>(py::str(dtype_py))};
  int num_elements = 1;
  for(int i = 0; i < ndim; i++) {
    num_elements *= a.get_shape()[i];
  }

	if (dtype == "float32") {
		return dot_impl<float>(
			num_elements, a.get_data(), b.get_data());
	} else if (dtype == "float64") {
		return dot_impl<double>(
			num_elements, a.get_data(), b.get_data());
	} else if (dtype == "int32") {
		return dot_impl<std::int32_t>(
			num_elements, a.get_data(), b.get_data());
	} else if (dtype == "int64") {
		return dot_impl<std::int64_t>(
			num_elements, a.get_data(), b.get_data());
	}
}

template<class T>
void print_flat_impl(int size, void * data){
	std::cout << "And your data values are: " << std::endl;
	for(int i = 0 ; i < size; i++)
		std::cout << reinterpret_cast<T*>(data)[i] << ", ";
	std::cout << std::endl;
}

void print_flat(const np::ndarray& a) {
  int ndim = py::extract<int>(a.attr("ndim"));
  auto shape = a.get_shape()[0];
  std::string dtype{py::extract<const char *>(py::str(a.get_dtype()))};
  auto num_elements = 1;
  for(int i = 0; i < ndim; i++) {
    num_elements *= a.get_shape()[i];
    std::cout << "Size of dim(" << i << "): " << a.get_shape()[i] << std::endl;
    std::cout << "Size of stride(" << i << "): " << a.get_strides()[i] << std::endl;
  }

	std::cout << "The type is: " << dtype << std::endl;

	if (dtype == "float32") {
		print_flat_impl<float>(num_elements, a.get_data());
	} else if (dtype == "float64") {
		print_flat_impl<double>(num_elements, a.get_data());
	} else if (dtype == "int32") {
		print_flat_impl<std::int32_t>(num_elements, a.get_data());
	} else if (dtype == "int64") {
		print_flat_impl<std::int64_t>(num_elements, a.get_data());
	}
}

BOOST_PYTHON_MODULE(cryptmat) {
  Py_Initialize();
  np::initialize();
  
  py::def("greet", greet);
  py::def("dot", dot);  
  py::def("print_flat", print_flat);
}
