#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>

#include "cryptly/test.hxx"

double L2(unsigned m, double x) { return 0.5*(x*x-2*(m+2)*x+(m+1)*(m+2)); }

constexpr auto wf = [](std::uint64_t L, std::uint64_t N, double x){
  return std::assoc_laguerre((2*L + 1), N - L - 1, double x) ;
};

TEST(Test, Hydrogen) {
  EXPECT_EQ(std::assoc_laguerre(2, 10, 0.5), L2(10, 0.5));
}
