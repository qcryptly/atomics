#include "cryptly/cuda_ops.hxx"
#include "cryptly/test.hxx"

TEST(Test, PrintDevices) { EXPECT_EQ(cryptly::print_devices(), 0); }
