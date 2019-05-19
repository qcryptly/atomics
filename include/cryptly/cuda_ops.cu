#include <iostream>
#include "cryptly/cuda_ops.hxx"

namespace cryptly {
inline namespace v1 {

int _r(int err) {
	if (err != cudaSuccess) {
		std::cout << "Error: " << cudaGetErrorString(cudaError_t(err)) << std::endl;
	}
	return err;
}

int print_devices() {
	int n_devices{};
	cudaError_t err = cudaGetDeviceCount(&n_devices);
	if (err != cudaSuccess) return _r(err);

	std::cout << "Number of devices we have: " << n_devices << std::endl;
	for (int i = 0; i < n_devices; i++) {
		cudaDeviceProp prop{};
		cudaGetDeviceProperties(&prop, i);
		std::cout << "Device name: " << prop.name << std::endl;
		std::cout << "Memory clock rate (KHz): " << prop.memoryClockRate << std::endl;
		std::cout << "Memory Bus Width (bits): " << prop.memoryBusWidth << std::endl;
		std::cout << "Peak Memory Bandwidth (GB/s): " << peak_mem_band(
			prop.memoryBusWidth, prop.memoryClockRate);
		std::cout << "Max Grid Size" << prop.maxGridSize << std::endl; 
		std::cout << "Max Threads Per Block" << prop.maxThreadsPerBlock << std::endl; 
		std::cout << "Max Threads Per Dim" << prop.maxThreadsDim << std::endl; 
		std::cout << "Warp Size" << prop.warpSize << std::endl; 
	}
	return cudaSuccess;
};

}
}
