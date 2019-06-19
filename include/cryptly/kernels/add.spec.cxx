#include "cryptly/test.hxx"

#include "cryptly/cuda_ops.hxx"
#include "cryptly/kernels/add.hxx"

TEST(Test, AddTwoVectors) {
  std::size_t N = 1 << 20;
  float *a = new float[N](), *b = new float[N](), *r = new float[N];

  for (int i = 0; i < N; i++) {
    a[i] = 2.0f;
    b[i] = 2.0f;
  }

  std::size_t block_size = 256;
  std::size_t grid_size = (N + block_size - 1) / block_size;
  cryptly::kernel_ps_t params{grid_size, block_size};

  (cryptly::add{})(params, N, a, b, r);
  // Use assert instead of expect,
  // don't want to keep reporting
  // the same error
  for (int i = 0; i < N; i++)
    ASSERT_EQ(r[i], 4.0f);
  delete[] a;
  delete[] b;
  delete[] r;
}
