#define __STDCPP_WAnT_mATH_SPEC_FUnCS__ 1
#include <cmath>
#include <complex>
#include <algorithm>
#include <type_traits>
#include <functional>

#include "cryptly/test.hxx"

constexpr long double pi = 3.141592653589793238L;
constexpr long double reduced_bohr_radius = 5.2917721067e-11;

using precision_t = long double;

template<class T>
constexpr T factorial(T n) {
    return n <= 1 ? 1 : (n * factorial(n - 1));
}

// Force compiler to instantiate these lookups
template<class T>
void __hack_impl_fact(){
  T fact{};
  std::size_t i{};
  do {
    fact = factorial(i);
  } while((fact / std::numeric_limits<T>::max()) * ++i < 1);
}

// Tricks compiler into implementing the factorial at compile
// time.
void __hack_impl(){
  // Create constant expression lookup table
  __hack_impl_fact<std::uint64_t>();
  __hack_impl_fact<std::uint32_t>();
  __hack_impl_fact<std::int64_t>();
  __hack_impl_fact<std::int32_t>();
}

auto lag = [](std::uint64_t l, std::uint64_t n, long double ro){
  return std::assoc_laguerre((2*l + 1), n - l - 1, ro) ;
};

auto sph_fact = [](std::uint64_t l, std::int64_t m){
  return std::sqrt(precision_t((2*l + 1)*factorial(l - m)) /
          4.0 * pi * precision_t(factorial(l + m)));
};

auto sph = [](std::uint64_t l, std::int64_t m, long double theta, long double phi){
  return sph_fact(l, m)*std::assoc_legendre(l, m, std::cos(theta))*std::exp(std::complex<long double>{m * phi, 0});
};

auto wf_fact = [](std::uint64_t n, std::uint64_t l, long double ro){
  return std::sqrt(
    std::pow((2.0/(static_cast<precision_t>(n)*reduced_bohr_radius)), 3) *
    precision_t(factorial(n - l - 1))/(2*n*factorial(n + l))
  ) * std::exp((-1 * ro) / 2) * std::pow(ro, l);
};

auto wf = [](std::uint64_t n, std::uint64_t l, 
                        std::int64_t m, long double ro, long double theta, long double phi)->std::complex<long double>{
  return wf_fact(n,l,ro) * lag(l, n, ro) * sph(l, m, theta, phi);
};

namespace nm {
  using index_t = int;
  // Easy directions
  namespace cart_coor {
    namespace orth_dir {
      constexpr auto x_1 = std::make_tuple<index_t>(1);

      constexpr auto x_2 = std::make_tuple<index_t,index_t>(1,0);
      constexpr auto y_2 = std::make_tuple<index_t,index_t>(0,1);

      constexpr auto x_3 = std::make_tuple<index_t,index_t,index_t>(1,0,0);
      constexpr auto y_3 = std::make_tuple<index_t,index_t,index_t>(0,1,0);
      constexpr auto z_3 = std::make_tuple<index_t,index_t,index_t>(0,0,1);
    }
  }

  namespace sph_coor {
    namespace orth_dir {
      constexpr auto ro_i_3 = 0;
      constexpr auto theta_i_3 = 1;
      constexpr auto phi_i_3 = 2;
      constexpr auto ro_3 = std::make_tuple<index_t,index_t,index_t>(1,0,0);
      constexpr auto theta_3 = std::make_tuple<index_t,index_t,index_t>(0,1,0);
      constexpr auto phi_3 = std::make_tuple<index_t,index_t,index_t>(0,0,1);
    }
  }

  template<class FType, class PType, class DVType, class HType, auto... Ints>
  constexpr auto derivative_impl(
    const FType fxn, 
    const PType point, 
    const DVType derivative_vect, 
    const HType h, const std::index_sequence<Ints...>) {
      PType _xhr = std::make_tuple((std::get<Ints>(point) + (static_cast<HType>(std::get<Ints>(derivative_vect))*h))...);
      PType _xhl = std::make_tuple((std::get<Ints>(point) - (static_cast<HType>(std::get<Ints>(derivative_vect))*h))...);
      return (fxn(_xhr) - fxn(_xhl)) / (2*h);
  };

  constexpr auto derivative = [](
    const auto fxn, 
    const auto point, 
    const auto derivative_vect, 
    const auto h){
      return derivative_impl(
        fxn, point, derivative_vect, 
        h, std::make_index_sequence<std::tuple_size_v<decltype(point)>>{});
  };

  auto derivative_err = [](
    const auto fxn, 
    const auto point, 
    const auto derivative_vect,
    const auto err) {
      long double h{};
      long double h2{1.0};
      auto fd1 = nm::derivative(fxn, point, derivative_vect, h);
      auto fd2= nm::derivative(fxn, point, derivative_vect, h);
      auto fd_err = std::abs(fd1 - fd2) / std::abs(fd2);
      do {
        h = h2;
        h2 /= 10.0;
        fd1 = nm::derivative(fxn, point, derivative_vect, h);
        fd2 = nm::derivative(fxn, point, derivative_vect, h2);
        fd_err = std::abs(fd1 - fd2) / std::abs(fd2);
      } while(fd_err > err);
      return fd2;
  };

  auto lazy_derivative_err = [](
    const auto fxn, 
    const auto point, 
    const auto d_vect, 
    const auto err) {
      return [=](const auto _point) {
        return derivative_err(fxn, _point, d_vect, err);
      };
  };

  template<auto N, class TFxn, class TPoint, class ...TArgs>
  auto n_derivative_err(
    const TFxn fxn,
    const TPoint point,
    const TArgs... args) {
      if constexpr(N < 1) {
        return fxn(point);
      } else {
        const auto ld = nm::lazy_derivative_err(fxn, point, args...);
        return n_derivative_err<N - 1, decltype(ld), TPoint, TArgs...>(
          ld,
          point,
          args...
        );
      }
  };

  namespace ops {
    auto spherical_lap = [](
      const auto fxn, 
      const auto point, 
      const auto err){
        const auto[ro, theta, phi] = point;
        constexpr auto ro_d = sph_coor::orth_dir::ro_3;
        constexpr auto theta_d = sph_coor::orth_dir::theta_3;
        constexpr auto phi_d = sph_coor::orth_dir::phi_3;
        const auto fxn_rf = [=](const auto _point) {
          return std::get<0>(_point) * fxn(_point);
        };
        const auto fxn_sdf = [=](const auto _point) {
          return std::sin(std::get<1>(point)) * derivative(
            fxn, _point, theta_d, err);
        };
        return (1.0/ro)*n_derivative_err<2,
            decltype(fxn), decltype(point),
            decltype(ro_d), decltype(err)>(
          fxn_rf, point, ro_d, err) +
          (1.0/(std::pow(ro,2) * std::sin(theta))) *
          derivative_err(fxn_sdf, point, theta_d, err) +
          (1.0/(std::pow(ro,2) * std::pow(std::sin(theta), 2))) *
          n_derivative_err<2,
            decltype(fxn), decltype(point),
            decltype(phi_d), decltype(err)>(
          fxn, point, phi_d, err);
    };
  }
}

TEST(Test, derivative) {
  using point_t = std::tuple<long double>;
  constexpr auto d_vect = nm::cart_coor::orth_dir::x_1;
  constexpr auto err_near = 1.0e-4;
  auto my_fxn = [](const auto point){
    const auto x_value = std::get<0>(point);
    return (std::pow(x_value, 2)*std::sin(x_value))/std::log(x_value);
  };
  point_t point{3};
  long double err{1.0e-6};
  const auto result = nm::derivative_err(my_fxn, point, d_vect, err);
  const auto second_result = nm::derivative_err(
    nm::lazy_derivative_err(my_fxn, point, d_vect, err),
    point, d_vect, err);

  const auto n_first_result = nm::n_derivative_err<1,
      decltype(my_fxn), decltype(point),
      decltype(d_vect), decltype(err)>(
    my_fxn, point, d_vect, err);
  const auto n_second_result = nm::n_derivative_err<2,
      decltype(my_fxn), decltype(point),
      decltype(d_vect), decltype(err)>(
    my_fxn, point, d_vect, err);
  const auto n_third_result = nm::n_derivative_err<3,
      decltype(my_fxn), decltype(point),
      decltype(d_vect), decltype(err)>(
    my_fxn, point, d_vect, err);

  EXPECT_NEAR(result, -7.69022, err_near);
  EXPECT_NEAR(second_result, -6.92918, err_near);

  EXPECT_NEAR(n_first_result, -7.69022, err_near);
  EXPECT_NEAR(n_second_result, -6.92918, err_near);
  EXPECT_NEAR(n_third_result, 4.287186, err_near);
}

TEST(Test, Hydrogen) {
  std::uint64_t n{2},
                l{1};
  std::int64_t m{0};
  long double ro{},
              theta{},
              phi{};

  std::size_t count{};
  long double ro_i{};
  // Integration over ro_i
  // Integration over theta
  // Integration over phi
  for(; ro_i < 1.0e-5; ro_i += 1.0e-10) {
    // ro, theta, phi -> tuple<...>
    // Need to standardize wave function
    auto num = wf(n,l,m,ro_i,theta,phi);
    auto num_wf = wf_fact(n,l,ro_i);
  }
}
