#pragma once

#include <complex>
#include <algorithm>
#include <type_traits>
#include <functional>

namespace cryptly::math {
inline namespace v1 {

using precision_t = long double;

// Throw the scientific constants up here
namespace physical_units {
    constexpr long double pi = 3.141592653589793238L;
    constexpr long double reduced_bohr_radius = 5.2917721067e-11;
    constexpr long double planks_constant_joules = 4.135667696e-15;
    constexpr long double reduced_planks_constant_joules = 1.054571817e-34; 
    constexpr long double electron_mass = 9.10938356e-31; // In kg 
    constexpr long double proton_mass = 1.6726219e-27; // In kg
    constexpr long double electron_charge = 1.602176634e-19; // In columbs
    constexpr long double free_perm = 8.8541878128e-12; // C^2 / N m2 --> N = kg m / s^2
    constexpr long double hydrogen_reduced_mass = (electron_mass * proton_mass) / (electron_mass + proton_mass);
    constexpr long double hamiltonian_potential_const = (electron_charge * electron_charge) / (4 * pi * free_perm);
    constexpr long double hamiltonian_kinetic_const = (reduced_planks_constant_joules * reduced_planks_constant_joules) / (2.0 * hydrogen_reduced_mass); 
}
////////////////////////////////////////////////////
// Define some common functions
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
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// This will work as our numerics lib...
// which includes useful directions and ops
namespace nm {
  using index_t = int;
  template <class T>
    struct is_complex {
	  typedef char one;
	  typedef struct { char x[2]; } two;
	  template <class C>
	  static one test (int C::*imag);
	  template <class C>
	  static two test (...);
	  public:
		enum {value = sizeof(is_complex<T>::test<T>(0)) == 1};
  };

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
      const auto diff = fxn(_xhr) - fxn(_xhl);
      if constexpr(is_complex<decltype(diff)>::value) {
        return diff / std::complex<decltype(diff.imag())>{2*h};
      } else {
        return diff / (2*h);
      }
  };

  template<class FType, class PType, class DVType, class HType, auto... Ints>
  constexpr auto derivative_forward_impl(
    const FType fxn, 
    const PType point, 
    const DVType derivative_vect, 
    const HType h, const std::index_sequence<Ints...>) {
      PType _xhr = std::make_tuple((std::get<Ints>(point) + (static_cast<HType>(std::get<Ints>(derivative_vect))*h))...);
      const auto diff = fxn(_xhr) - fxn(point);
      if constexpr(is_complex<decltype(diff)>::value) {
        return diff / std::complex<decltype(diff.imag())>{h};
      } else {
        return diff / h;
      }
  };

  template<bool Central = true>
  constexpr auto derivative = [](
    const auto fxn, 
    const auto point, 
    const auto derivative_vect, 
    const auto h){
      if constexpr(Central) {
          return derivative_impl(
            fxn, point, derivative_vect, 
            h, std::make_index_sequence<std::tuple_size_v<decltype(point)>>{});
      } else {
          return derivative_forward_impl(
            fxn, point, derivative_vect, 
            h, std::make_index_sequence<std::tuple_size_v<decltype(point)>>{});
      }
  };

  // We need a concept for complex conjugate vs real numbers
  template<bool Central = true>
  auto derivative_err = [](
    const auto fxn, 
    const auto point, 
    const auto derivative_vect,
    const auto err) {
      long double h{};
      long double h2{1.0};
      auto fd1 = nm::derivative<Central>(fxn, point, derivative_vect, h);
      auto fd2= nm::derivative<Central>(fxn, point, derivative_vect, h);
      auto fd_err = std::abs(fd1 - fd2) / std::abs(fd2);
      do {
        h = h2;
        h2 /= 10.0;
        fd1 = nm::derivative<Central>(fxn, point, derivative_vect, h);
        fd2 = nm::derivative<Central>(fxn, point, derivative_vect, h2);
        fd_err = std::abs(fd1 - fd2) / std::abs(fd2);
      } while(fd_err > err);
      return fd2;
  };

  template<bool Central = true>
  auto lazy_derivative_err = [](
    const auto fxn, 
    const auto point, 
    const auto d_vect, 
    const auto err) {
      return [=](const auto _point) {
        return derivative_err<Central>(fxn, _point, d_vect, err);
      };
  };

  template<bool Central = true, auto N, class TFxn, class TPoint, class ...TArgs>
  auto n_derivative_err(
    const TFxn fxn,
    const TPoint point,
    const TArgs... args) {
      if constexpr(N < 1) {
        return fxn(point);
      } else {
        const auto ld = nm::lazy_derivative_err<Central>(fxn, point, args...);
        return n_derivative_err<Central, N - 1, decltype(ld), TPoint, TArgs...>(
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
        const auto fxn_rf = [_fxn = fxn](const auto _point) {
          return std::get<sph_coor::orth_dir::ro_i_3>(_point) * _fxn(_point);
        };
        const auto fxn_sdf = [_fxn = fxn, _theta_d = theta_d, _err = err](const auto _point) {
          return std::sin(std::get<sph_coor::orth_dir::theta_i_3>(_point)) * derivative<false>(
            _fxn, _point, _theta_d, _err);
        };
        return (1.0/ro)*n_derivative_err<false, 2,
            decltype(fxn_rf), decltype(point),
            decltype(ro_d), decltype(err)>(
          fxn_rf, point, ro_d, err) +
          (1.0/(std::pow(ro,2) * std::sin(theta))) *
          derivative_err<true>(fxn_sdf, point, theta_d, err) +
          (1.0/(std::pow(ro,2) * std::pow(std::sin(theta), 2))) *
          n_derivative_err<false, 2,
            decltype(fxn), decltype(point),
            decltype(phi_d), decltype(err)>(
          fxn, point, phi_d, err);
    };

    auto hamiltonian = [](
      const auto wavefunction, 
      const auto point, 
      const auto err){
        const auto[ro, theta, phi] = point;
        constexpr auto ro_d = sph_coor::orth_dir::ro_3;
        constexpr auto theta_d = sph_coor::orth_dir::theta_3;
        constexpr auto phi_d = sph_coor::orth_dir::phi_3;
        return (-1)*physical_units::hamiltonian_kinetic_const * spherical_lap(wavefunction, point, err) - physical_units::hamiltonian_potential_const * (1.0/ro) * wavefunction(point);
    };
  }
}
// End of ops namespace
namespace models {
	using scalar_t = long double;
	using point_t = std::tuple<scalar_t, scalar_t, scalar_t>;
	auto lag = [](std::uint64_t n, std::uint64_t l, long double ro){
	  return std::assoc_laguerre((2*l + 1), n - l - 1, ro) ;
	};

	auto sph_fact = [](std::uint64_t l, std::int64_t m){
	  return std::sqrt(precision_t((2*l + 1)*factorial(l - m)) /
			  4.0 * physical_units::pi * precision_t(factorial(l + m)));
	};

	auto sph = [](std::uint64_t l, std::int64_t m, long double theta, long double phi){
	  return sph_fact(l, m)*std::assoc_legendre(l, m, std::cos(theta))*std::exp(std::complex<long double>{m * phi, 0});
	};

	auto wf_fact = [](std::uint64_t n, std::uint64_t l, long double ro){
	  return std::sqrt(
		std::pow((2.0/(static_cast<precision_t>(n)*physical_units::reduced_bohr_radius)), 3) *
		precision_t(factorial(n - l - 1))/(2*n*factorial(n + l))
	  ) * std::exp((-1 * ro) / 2) * std::pow(ro, l);
	};

    // Simple electron model wave function
    // W(ro) * W(theta, phi)
	auto wf_sem_simple = [](const auto hyper_param) {
		const auto[_n,_l,_m] = hyper_param;
		return [n = _n, l = _l, m = _m](point_t point){
			const auto[ro, theta, phi] = point;	
			return wf_fact(n,l,ro) * lag(n, l, ro) * sph(l, m, theta, phi);
		};
	};
}

// End of models namespace

}
}
