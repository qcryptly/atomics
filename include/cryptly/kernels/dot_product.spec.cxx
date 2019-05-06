#include "cyrptly/test.hxx"
#include "cryptly/kernels/dot_product.hxx"

TEST(dot_product, DotProductSync) {
	constexpr std::uint64_t size_test = 1 << 20;
	constexpr auto block_size = 256;

	// This is a weird value to set.
	// If it's too large you will seg fault
	// Must know what the device specs are.
	// Lets just set to 256 to keep things simple
	constexpr auto grid_size = 256;
	kernel_ps_t params{block_size, grid_size};	
	kernel_ps_t kernel{block_size, grid_size};


	// Allocate on Host
	EXPECT_EQ(product, (size_test)*4);
	return 0;
}
