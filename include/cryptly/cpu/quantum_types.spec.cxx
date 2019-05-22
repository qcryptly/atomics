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
