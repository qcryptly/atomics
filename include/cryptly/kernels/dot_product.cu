#include "cryptly/cuda_ops.hxx"
#include "cryptly/kernels/dot_product.hxx"

template<int BSize = 256, class TArg>
__global__
void dot_product_cuda(std::size_t N, const TArg* a, const TArg* b, TArg* r){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	__shared__ TArg agg[BSize];
	agg[threadIdx.x] = 0;
	// Map
	for (int i = index; i < N; i += stride) {
		agg[threadIdx.x] += a[i] * b[i];
	}

	// Wait for all block threads to finish
	__syncthreads();

	// Reduce everything
	// Use one writer,
	// warp divergence
	if (threadIdx.x == 0) {
		r[blockIdx.x] = 0; 
		for (int i = 0; i < blockDim.x; i ++) {
			r[blockIdx.x] += agg[i];
		}
	}
}

template<class T, typename... Ts>
int dot_product_impl<T, Ts...>::operator()(const kernel_ps_t& ps,
	std::size_t N, const T* a_i, const T* b_i, T* r_o) {
	T *cuda_mem_a, *cuda_mem_b, *cuda_mem_r;
	auto size_n = sizeof(T) * N;
	auto size_r = sizeof(T) * ps.grid_size;
	// Allocate on Device
	_r(cudaMalloc(&cuda_mem_a, size_n));	
	_r(cudaMalloc(&cuda_mem_b, size_n));	
	_r(cudaMalloc(&cuda_mem_r, size_r));	

	// Copy Host to Device
	_r(cudaMemcpy(cuda_mem_a, a_i, size_n, cudaMemcpyHostToDevice));
	_r(cudaMemcpy(cuda_mem_b, b_i, size_n, cudaMemcpyHostToDevice));

	// Run Kernel
	// Hard code grid size for now
	dot_product_cuda<256, T><<<ps.grid_size, ps.block_size, size_r, (CUstream_st*)ps.cuda_stream>>>(N, cuda_mem_a, cuda_mem_b, cuda_mem_r);	
	
	cudaDeviceSynchronize();
	// Copy Device to Host
	T results[ps.grid_size];
	_r(cudaMemcpy(results, cuda_mem_r, size_r, cudaMemcpyDeviceToHost));

	// Simple aggregation
	T product{};
	for(int i = 0; i < ps.grid_size; i++) {
		product += results[i];
	}
	// Serialize and print results
	cudaFree(cuda_mem_a);
	cudaFree(cuda_mem_b);
	cudaFree(cuda_mem_r);
	*r_o = product;
	return cudaGetLastError();
}

// This is a hack to force the compiler
// to instantiate add_impl. If you leave
// this out, you will get undefined symbol
// reference for add_impl member functions :(
template <typename TArg, typename... Ts>
void __hack_impl_dot_product__(){
	[[maybe_unused]]TArg* _ = nullptr;
	(dot_product{})(
		kernel_ps_t{}, std::size_t{},
		_, _, _); 
	__hack_impl_dot_product__<Ts...>();
}

template<>
void __hack_impl_dot_product__<bool>(){
	[[maybe_unused]]bool* _ = nullptr;
	(dot_product{})(
		kernel_ps_t{}, std::size_t{},
		_, _, _);
}

void __hack_dot_product__(){
	__hack_impl_dot_product__<primitive_types>();
}
