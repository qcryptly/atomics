#include <boost/rational.hpp>
#include <cmath>
#include <iostream>
#include <sstream>
#include <tuple>

#include "cryptly/test.hxx"
#include "cryptly/cpu/quantum_types.hxx"

TEST(Test, MakeQuantumTypes) {
	cryptly::state_s state{{4,2},2,3,1};

	cryptly::jm_s jm_a{.j={1,2},.m={1,2}};
	cryptly::jm_s jm_b{.j={1,2},.m={-1,2}};
	// Assume 3 is J and m = 2
	cryptly::cg <double>table{jm_a,jm_b};
	EXPECT_EQ(table({.j={1,1},.m={0}}), std::sqrt(1.0/2.0));
}

TEST(QuantumTypes, CrossProduct) {
	using state_t = std::int64_t;
	using value_t = std::int64_t;
	using wf_vpair_t = std::pair<state_t, value_t>;
	using wf_vect_t = std::vector<wf_vpair_t>;
	wf_vect_t wf_a {
		{1,0}, {2,1}, {3,0}
	};
	wf_vect_t wf_b {
		{1,1}, {2,0}, {3,0}
	};
	cryptly::matrix_s<state_t, value_t> wf_mat{};
	wf_mat.push_back(wf_a);
	wf_mat.push_back(wf_b);

	wf_vect_t expected {
		{1,0}, {2,0}, {3,-1}
	};
	wf_vect_t actual {
		{1,0}, {2,0}, {3,0}
	};
	cryptly::cross<value_t>(wf_mat, actual);
	EXPECT_EQ(expected, actual);
}

TEST(QuantumTypes, Determinant) {
	using state_t = std::int64_t;
	using value_t = std::int64_t;
	using wf_vpair_t = std::pair<state_t, value_t>;
	using wf_vect_t = std::vector<wf_vpair_t>;
	wf_vect_t wf_a{
		{1,0}, {2,1}
	};
	wf_vect_t wf_b{
		{1,1}, {2,0}
	};
	value_t expected = -1;
	cryptly::matrix_s<state_t, value_t> wf_mat{};
	wf_mat.push_back(wf_a);
	wf_mat.push_back(wf_b);
	auto actual = cryptly::det<value_t>(wf_mat);
	EXPECT_EQ(actual, value_t{-1});
}

TEST(QuantumTypes, JDown) {
	std::stringstream buf{};
	cryptly::jm_s init{{5,2},{5,2}};
	auto starting_state = std::make_tuple(init, init);
	cryptly::wfj_t wfa {{}, cryptly::wfj_hash};
	wfa[starting_state] = cryptly::ev_vect_t{};
	cryptly::ev_vect_t init_ev {{}, cryptly::rational_hash};
	init_ev[1] = 1;
	wfa[starting_state] = init_ev;
	for (int limit = 5; limit > -1; limit --) {
		std::cout << "Printing M: " << limit << std::endl;
		cryptly::print_wf(wfa, buf);
		std::cout << buf.str();
		buf.str("");
		if (limit == 0) return;
		wfa = cryptly::j_down(wfa,
			cryptly::j_t {5,1},
			cryptly::m_t {limit,1});
	}
}

constexpr auto vtwm = [](const auto& v){
	cryptly::wfj_t wfa{{}, cryptly::wfj_hash};
	for (auto& e : v) wfa[e.first] = e.second;
	return wfa;
};

constexpr auto mtwv = [](const auto& v){
	using wfv = std::vector<std::pair<cryptly::wfj_tup_t, cryptly::ev_vect_t>>;
	return wfv{v.begin(), v.end()};
};

// We assume that we want to return m - 1
constexpr auto cg_m = [](
		const auto& initial_state,
		const auto& j,
		const auto& m
	){
	std::stringstream buf{};
	using matrix_t = typename std::remove_reference<decltype(initial_state)>::type;
	//using key_t = matrix_t::key_t;
	//using val_t = matrix_t::val_t;
	matrix_t spectrum{};
	if (m == j) {
		for (auto wf : initial_state._matrix) {
			auto first = cryptly::j_down(
				vtwm(wf),
				j,
				m
			);
			auto second = mtwv(first);
			// Handle the edge case that we are only going down the first level
			(second[0].second.begin())->second = -1 * (second[0].second.begin())->second;
			//spectrum.push_back(first);
			//spectrum.push_back(second);
			return spectrum;
		}
	}
	
	// Here we need to create a table that is from
	// j - m = size
	for(auto index = j; index > m - 1; index--) {
		// Calculate the the down states
		for (auto wf : initial_state._matrix) {
			/*spectrum.push_back(
				mtwv(cryptly::j_down(
					vtwm(wf), j - index, m
				))
			);*/
		}
//		spectrum.push_back(
//			cryptly::cross<cryptly::ev_vect_t>(wf_mat, actual));
	}
};

TEST(QuantumTypes, JDownCross) {
	std::stringstream buf{};
	
	// Initialize coupled state
	cryptly::jm_s init{{5,2},{5,2}};
	auto starting_state = std::make_tuple(init, init);

	// Create our wave function
	const cryptly::wfj_t wfa_con {{}, cryptly::wfj_hash};
	cryptly::wfj_t wfa {{}, cryptly::wfj_hash};
	
	// Build out our eigenvalue, root vector
	wfa[starting_state] = cryptly::ev_vect_t{};
	cryptly::ev_vect_t init_ev {{}, cryptly::rational_hash};
	init_ev[1] = 1; // (1)sqrt(0)
	wfa[starting_state] = init_ev;

	// Create our initial matrix
	using mpair_type = std::pair<cryptly::wfj_tup_t, cryptly::ev_vect_t>;
	using matrix_t = cryptly::matrix_s<cryptly::wfj_tup_t, cryptly::ev_vect_t>;
	matrix_t wf_mat{};
	wf_mat.push_back(matrix_t::state_t{wfa.begin(), wfa.end()});
	cg_m(wf_mat, cryptly::j_t{5,1}, cryptly::m_t{5,1});	
}
