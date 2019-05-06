#include "cryptly/test.hxx"

#include "cryptly/cuda_ops.hxx"
#include "cryptly/kernels/add.hxx"

TEST(add, AddTwoVectors) {
	int N = 1<<20;
	float *a, *b, *r;

	for (int i = 0; i < N; i++) {
		a[i] = 2.0f;
		b[i] = 2.0f;
	}

	// User entire query

	int block_size = 256;
	int grid_size = (N + block_size - 1) / block_size;

	kernel_ps_t params{grid_size, block_size};
	_r((add{})(params, N, a, b));

	EXPECT_EQ(maxError, 0.0f);
	return 0;
}
