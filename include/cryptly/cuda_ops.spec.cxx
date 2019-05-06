#include "cryptly/test.hxx"
#include "cryptly/cuda_ops.hxx"

TEST(Test, PrintDevices) {
	EXPECT_EQ(print_devices(), 0);
}
