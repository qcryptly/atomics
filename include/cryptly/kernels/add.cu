#include "cryptly/cuda_ops.hxx"
#include "cryptly/kernels/add.hxx"

// Kernel function to add the elements of two arrays
template<class TArg=int>
__global__
void add_cuda(std::size_t n, TArg *a, TArg *b, TArg *r)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride)
    r[i] = a[i] + b[i];
}

template<typename T, typename... Ts>
int add_impl<T, Ts...>::operator()(const kernel_ps_t& ps,
	std::size_t N, const T* a_i, const T* b_i, T* r){
	cudaError_t error{};
	T *cuda_mem_a, *cuda_mem_b, *cuda_mem_r;
	auto n_bytes = N * sizeof(T);
	// _r logs any cuda errors, can turn this off
	// in production
	_r(cudaMalloc(&cuda_mem_a, n_bytes));
	_r(cudaMalloc(&cuda_mem_b, n_bytes));
	_r(cudaMalloc(&cuda_mem_r, n_bytes));
	// Copy bytes from a_i and b_i
	_r(cudaMemcpy(cuda_mem_a, a_i, n_bytes, cudaMemcpyHostToDevice));	
	_r(cudaMemcpy(cuda_mem_b, b_i, n_bytes, cudaMemcpyHostToDevice));	
	add_cuda<T><<<ps.grid_size, ps.block_size, ps.shared_size, (CUstream_st*)ps.cuda_stream>>>(N, cuda_mem_a, cuda_mem_b, cuda_mem_r);
	cudaDeviceSynchronize();
	if (error == cudaSuccess)
		error = _r(cudaMemcpy(r, cuda_mem_r, n_bytes, cudaMemcpyDeviceToHost));	
	cudaFree(cuda_mem_a);
	cudaFree(cuda_mem_b);
	return error;	
}
