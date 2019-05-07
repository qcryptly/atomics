#include "cryptly/test.hxx"
#include "cryptly/kernels/dot_product.hxx"

TEST(dot_product, DotProductSync) {
	constexpr std::uint64_t N = 1 << 20;
	float *a = new float[N](),
		  *b = new float[N](),
		  r;
	for (int i = 0; i < N; i++)
		{a[i] = 2.0f; b[i] = 2.0f;}	

	constexpr auto block_size = 256;

	// This is a weird value to set.
	// If it's too large you will seg fault
	// Must know what the device specs are.
	// Lets just set to 256 to keep things simple
	constexpr auto grid_size = 256;
	kernel_ps_t kernel{block_size, grid_size};
	
	(dot_product{})(kernel, N, a, b, &r);

	// Allocate on Host
	ASSERT_EQ(r, (N)*4);
	delete[] a;
	delete[] b;
}
