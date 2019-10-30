#include <tuple>
#include <iostream>

#include "cryptly/math/models.hxx"

#include "cryptly/test.hxx"

TEST(Test, derivative) {
	using namespace cryptly::math;
	using point_t = std::tuple<long double>;
	constexpr auto d_vect = nm::cart_coor::orth_dir::x_1;
	constexpr auto err_near = 1.0e-4;
	constexpr auto err{1.0e-6};
	point_t point{3};

	auto my_fxn = [](const auto point){
		const auto x_value = std::get<0>(point);
		return (std::pow(x_value, 2)*std::sin(x_value))/std::log(x_value);
	};
	const auto result = nm::derivative<true>(my_fxn, point, d_vect, err); 
	const auto second_result = nm::derivative_err<true>(
		nm::lazy_derivative_err<true>(my_fxn, point, d_vect, err),
		point, d_vect, err);

	const auto n_first_result = nm::n_derivative_err<true, 1,
		decltype(my_fxn), decltype(point),
		decltype(d_vect), decltype(err)>(
			my_fxn, point, d_vect, err);
	const auto n_second_result = nm::n_derivative_err<true, 2,
		decltype(my_fxn), decltype(point),
		decltype(d_vect), decltype(err)>(
			my_fxn, point, d_vect, err);
	const auto n_third_result = nm::n_derivative_err<true, 3,
		decltype(my_fxn), decltype(point),
		decltype(d_vect), decltype(err)>(
			my_fxn, point, d_vect, err);

	EXPECT_NEAR(result, -7.69022, err_near);
	EXPECT_NEAR(second_result, -6.92918, err_near);

	EXPECT_NEAR(n_first_result, -7.69022, err_near);
	EXPECT_NEAR(n_second_result, -6.92918, err_near);
	EXPECT_NEAR(n_third_result, 4.287186, err_near);
}

TEST(Test, SingleElectronModel) {
	using namespace cryptly::math;
    auto ham_density_fxn = [](
            const auto quantum_numbers,
            const auto coordinates,
            const auto err
        ){
        const auto [ro, theta, phi] = coordinates;
        auto wf = models::wf_sem_simple(quantum_numbers);
        const auto result = wf(coordinates);
        const auto op_results = nm::ops::spherical_lap(wf, coordinates, err);
        const auto ham_results = nm::ops::hamiltonian(wf, coordinates, err);
        const auto expect = std::conj(wf(coordinates)) * ham_results;
        return expect * std::pow(ro, 2) * std::sin(phi);
    };
	constexpr auto err{1.0e-6};
	// n, l, m
	const auto quantum_numbers = std::make_tuple(2,1,0); 
	const auto coordinates = std::make_tuple<long double, long double, long double>(1.0, 0.4, 0.3);
    std::cout << "HDF: " << ham_density_fxn(quantum_numbers, coordinates, err) << std::endl;
//    const auto [ro, theta, phi] = coordinates;
//	auto wf = models::wf_sem_simple(quantum_numbers);
//	const auto result = wf(coordinates);
//	std::cout << "Wave function: " << result << std::endl;
//	const auto op_results = nm::ops::spherical_lap(wf, coordinates, err);
//	std::cout << "Spherical Lap applied: " << op_results << std::endl;
//    const auto ham_results = nm::ops::hamiltonian(wf, coordinates, err);
//    std::cout << "Hamiltonian applied: " << ham_results << std::endl;
//    const auto expect = std::conj(wf(coordinates)) * ham_results;
//    std::cout << "Expectation: " << expect << std::endl;
//    const auto expect_space = expect * std::pow(ro, 2) * std::sin(phi);
//    std::cout << "Expected space: " << expect_space << std::endl;
  // n, l, m, ro, theta, phi
//  using fxns = cryptly::math::models::fxns;
//  using ops = cryptly::math::models::ops;
//  const auto electron_numbers = std::make_tuple(2,1,0);
//  const auto electron_position = std::make_tuple(0.0,0.0,0.0);
//  // Lazy lambda
//  const auto sem = fxns::sem(qm_numbers);
//  // Conjugate sem
//  const auto sem_i = ops::conj(fxns::sem, electron_numbers);
//  const auto ops_sph_lap = ops::sph_lap;
//  long double err = 1.0e-6;
//  std::complex<long double> expected_value{0,0};
//  const auto res_1 = sph_lap(sem, electron_position, err);
//  // Due to eager computation, we can just use res_1;
//  const auto res_2 = sem_i(electron_position) * res_1;
//  EXPECT_EQ(sem(electron_position), expected_value);
//  EXPECT_EQ(res_1, expected_value); 
//  EXPECT_EQ(res_2, expected_value); 
}
