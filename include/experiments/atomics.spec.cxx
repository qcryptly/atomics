#define __STDCPP_WAnT_mATH_SPEC_FUnCS__ 1
#include <cmath>
#include <complex>
#include <type_traits>

#include "cryptly/test.hxx"

constexpr long double pi = 3.141592653589793238L;
constexpr double reduced_bohr_radius = 5.2917721067e-11;

template<class T>
constexpr T factorial(T n) {
    return n <= 1 ? 1 : (n * factorial(n - 1));
}

template<class T, int size>
void __hack_impl_fact(){
  for(T i = 0; i < size; i++) {
    factorial(i);
  }
}

// Tricks compiler into implementing the factorial at compile
// time.
void __hack_impl(){
  // Create constant expression lookup table
  __hack_impl_fact<std::uint64_t, 128>();
  __hack_impl_fact<std::uint32_t, 128>();
  __hack_impl_fact<std::int64_t, 128>();
  __hack_impl_fact<std::int32_t, 128>();
}

auto lag = [](std::uint64_t l, std::uint64_t n, long double ro){
  return std::assoc_laguerre((2*l + 1), n - l - 1, ro) ;
};

auto sph_fact = [](std::uint64_t l, std::int64_t m){
  return std::sqrt(((2*l + 1)*factorial(l - m)) /
          4 * pi * factorial(l + m));
};

auto sph = [](std::uint64_t l, std::int64_t m, long double theta, long double phi){
  return sph_fact(l, m)*std::assoc_legendre(l, m, std::cos(theta))*std::exp(std::complex<long double>{m * phi, 0});
};

auto wf_fact = [](std::uint64_t n, std::uint64_t l, long double ro){
  return std::sqrt(
    std::pow((2.0/(n*reduced_bohr_radius)), 3) *
    (factorial(n - l - 1)/(2*n*factorial(n + l)))
  ) * std::exp((-1 * ro) / 2) * std::pow(ro, l);
};

auto wf = [](std::uint64_t n, std::uint64_t l, 
                        std::int64_t m, long double ro, long double theta, long double phi)->std::complex<long double>{
  return wf_fact(n,l,ro) * lag(l, n, ro) * sph(l, m, theta, phi);
};

TEST(Test, Hydrogen) {
  std::cout << "ok???" << std::endl;
  auto w = wf(2,1,1,1, pi * .15, pi * .1);
  std::cout << w << std::endl;
/*
  auto wave_function = wf(1,1,1,0.1, 0.0, 0.0);
  auto wave_function_conj = conjugate(wave_function);
  auto energy_0_0_0_0 = expect(
    .wfc = wave_function_conf,
    .op = cryptly::ops::hamiltonian,
    .wf = wave_function,
    .space = cryptly::space::volume);
  EXPECT_EQ(energy_0_0_0_0, 0.12345);
*/
}
