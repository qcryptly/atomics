#include <iostream>

#include "cryptly/test.hxx"

#include "cryptly/math/rational.hxx"
#include "cryptly/math/linear_ops.hxx"

using number_t = std::int64_t;
using component_t = typename cryptly::math::rational_s<number_t>::component_t;
using hashmap_t = typename cryptly::math::rational_s<number_t>::hashmap_t;
using hashkey_t = typename cryptly::math::rational_s<number_t>::hashkey_t;
using hashval_t = typename cryptly::math::rational_s<number_t>::hashval_t;
using vector_t = typename cryptly::math::rational_s<number_t>::vector_t;
using hash_t = typename cryptly::math::rational_s<number_t>::hash_t;
using base_t = typename cryptly::math::rational_s<number_t>::base_t;
using rational_t = typename cryptly::math::rational_s<number_t>; 

using matrix_n_t = typename cryptly::math::matrix_s<number_t, number_t>;
using _matrix_n_t = typename matrix_n_t::matrix_t;
using state_t = typename matrix_n_t::state_t;

using matrix_d_t = typename cryptly::math::matrix_s<double, double>;
using _matrix_d_t = typename matrix_d_t::matrix_t;
using state_d_t = typename matrix_d_t::state_t;

using stater_t = typename cryptly::math::matrix_s<rational_t, rational_t>::state_t;

constexpr double coarse_eps = 1.0e-3;
constexpr double bad_eps = 1.0e-6;
constexpr double low_eps = 1.0e-9;
constexpr double med_eps= 1.0e-14;
constexpr double high_eps= std::numeric_limits<double>::epsilon();

TEST(LinearOps, Cross) {
	// Normal numbers
	cryptly::math::matrix_s<number_t, number_t> spectrum{};
	state_t a{{0,0},{1,1},{2,0}};
	state_t b{{0,1},{1,0},{2,0}};
	state_t expected{{0,0},{1,0},{2,-1}};	
	state_t actual{{0,0},{1,0},{2,0}};	
	spectrum.push_back(a);
	spectrum.push_back(b);
	spectrum.cross(actual);

	EXPECT_EQ(expected, actual);

	// Crytply rationals	
	cryptly::math::matrix_s<rational_t, rational_t> spectrumr{};
	hashmap_t zero_hash {{
		{{1,1}, 0}
	}, 1, cryptly::math::hashkey_rational_h<hashkey_t>};
	rational_t zero{zero_hash}; 
	hashmap_t one_hash {{
		{{1,1}, 1}
	}, 1, cryptly::math::hashkey_rational_h<hashkey_t>};
	rational_t one{one_hash};
	hashmap_t two_hash {{
		{{1,1}, 2}
	}, 1, cryptly::math::hashkey_rational_h<hashkey_t>};
	rational_t two{two_hash};
	auto inverse = one * -1;
	stater_t ar{{zero, zero}, {one, one}, {two, zero}};
	stater_t br{{zero, one}, {one, zero}, {two, zero}};
	stater_t expectedr{{zero, zero}, {one, zero}, {two, inverse}};
	stater_t actualr{{zero, zero}, {one, zero}, {two, zero}};
	spectrumr.push_back(ar);	
	spectrumr.push_back(br);	
	spectrumr.cross(actualr);	
	EXPECT_EQ(expectedr, actualr);
}

TEST(LinearOps, Orth) {
	matrix_d_t spectrum_a{{
		state_d_t{{0,1},{1,3},{2,0}},
		state_d_t{{0,2},{1,6},{2,5}},
		state_d_t{{0,3},{1,9},{2,1}},
	}};
	matrix_d_t expected_a{{
		state_d_t{{0,1},{1,3},{2,0}},
		state_d_t{{0,0},{1,0},{2,1}},
		state_d_t{{0,0},{1,0},{2,0}},
	}};
	auto actual_a = matrix_d_t{spectrum_a.rref()};
	EXPECT_EQ(actual_a, expected_a);

	matrix_d_t spectrum_b{{
		state_d_t{{0,0},{1,3},{2,-6},{3,6},{4,4},{5,-5}},
		state_d_t{{0,3},{1,-7},{2,8},{3,-5},{4,8},{5,9}},
		state_d_t{{0,3},{1,-9},{2,12},{3,-9},{4,6},{5,15}},
	}};
	spectrum_b.eps = med_eps;
	matrix_d_t expected_b{{
		state_d_t{{0,1},{1,0},{2,-2},{3,3},{4,0},{5,-24}},
		state_d_t{{0,0},{1,1},{2,-2},{3,2},{4,0},{5,-7}},
		state_d_t{{0,0},{1,0},{2,0},{3,0},{4,1},{5,4}},
	}};
	auto actual_b = matrix_d_t{spectrum_b.rref()};
	actual_b.eps = low_eps;
	EXPECT_EQ(actual_b, expected_b);

	matrix_d_t spectrum_c{{
		state_d_t{{0,0},{1,0},{2,5},{3,-8},{4,2},{5,47}},
		state_d_t{{0,0},{1,0},{2,6},{3,4},{4,6},{5,-7}},
		state_d_t{{0,0},{1,0},{2,0},{3,0},{4,0},{5,1}},
		state_d_t{{0,2},{1,5},{2,6},{3,-4},{4,7},{5,4}},
	}};
	spectrum_c.eps = med_eps;
	matrix_d_t expected_c{{
		state_d_t{{0,1},{1,2.5},{2,0},{3,0},{4,1.55882},{5,0}},
		state_d_t{{0,0},{1,0},{2,1},{3,0},{4,0.823529},{5,0}},
		state_d_t{{0,0},{1,0},{2,0},{3,1},{4,0.264706},{5,0}},
		state_d_t{{0,0},{1,0},{2,0},{3,0},{4,0},{5,1}},
	}};
	
	auto actual_c = matrix_d_t{spectrum_c.rref()};
	actual_c.eps = coarse_eps;
	EXPECT_EQ(actual_c, expected_c);

	// cg
	matrix_d_t spectrum_cg_init{{
		state_d_t{{0,0.707107},{1,0.707107}},
	}};
	spectrum_cg_init.eps = high_eps;

	matrix_d_t expected_cg_init{{
		state_d_t{{0,0.707107},{1,-0.707107}},
	}};
	for(auto x : spectrum_cg_init.cross()) {
		std::cout << x.second << ", " ;
	}
	
//	auto actual_cg = matrix_d_t{{spectrum_cg_init.cross()}};
//	actual_cg.eps = high_eps;
//	EXPECT_EQ(actual_cg, expected_cg);

	// cg
	matrix_d_t spectrum_cg{{
		state_d_t{{0,0.288675135},{1,0.645497224},{2,0.645497224},{3,0.288675135}},
		state_d_t{{0,0.56694671},{1,0.422577127},{2,-0.422577127},{3,-0.56694671}},
		state_d_t{{0,0.645497224},{1,-0.288675135},{2,-0.288675135},{3,0.645497224}},
	}};

	matrix_d_t expected_cg{{
		state_d_t{{0,0.422577},{1,-0.566947},{2,0.566947},{3,-0.422577}},
	}};
	
	auto actual_cg = matrix_d_t{{spectrum_cg.cross()}};
	actual_cg.eps = bad_eps;
	EXPECT_EQ(actual_cg, expected_cg);
}
