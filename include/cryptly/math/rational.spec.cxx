#include <iostream>
#include "cryptly/test.hxx"
#include "cryptly/math/rational.hxx"

using component_t = typename cryptly::math::rational_s<std::int64_t>::component_t;
using hashmap_t = typename cryptly::math::rational_s<std::int64_t>::hashmap_t;
using hashkey_t = typename cryptly::math::rational_s<std::int64_t>::hashkey_t;
using hashval_t = typename cryptly::math::rational_s<std::int64_t>::hashval_t;
using vector_t = typename cryptly::math::rational_s<std::int64_t>::vector_t;
using hash_t = typename cryptly::math::rational_s<std::int64_t>::hash_t;
using base_t = typename cryptly::math::rational_s<std::int64_t>::base_t;
	
TEST(Rational, EigenRootsAddition) {
	hashmap_t initial_hash{{
		{{4,{1,2}}, 2},
		{{16,{1,2}}, 2}
	}, 2, cryptly::math::hashkey_rational_h<hashkey_t>};

	hashmap_t expected_hash{{
		{{4,{1,2}}, 4},
		{{16,{1,2}}, 4}
	}, 2, cryptly::math::hashkey_rational_h<hashkey_t>};

	cryptly::math::rational_s<std::int64_t> roots_a{initial_hash};
	cryptly::math::rational_s<std::int64_t> roots_b{initial_hash};
	cryptly::math::rational_s<std::int64_t> expected{expected_hash};
	
	auto roots_c = roots_a + roots_b;
	EXPECT_EQ(roots_a, roots_b);
	EXPECT_NE(roots_a, roots_c);
	EXPECT_EQ(expected, roots_c);
}

TEST(Rational, EigenRootsMultiplication) {
	hashmap_t initial_hash{{
		{{4,{1,2}}, 2},
		{{16,{1,2}}, 2}
	}, 2, cryptly::math::hashkey_rational_h<hashkey_t>};

	cryptly::math::rational_s<std::int64_t> roots_a{initial_hash};
	cryptly::math::rational_s<std::int64_t> roots_b{initial_hash};
	auto roots_c = roots_a * roots_b;
	EXPECT_EQ(roots_a, roots_b);
	EXPECT_NE(roots_a, roots_c);
	EXPECT_EQ(double(roots_c), 144);
	roots_c = roots_c * 2;
	EXPECT_EQ(double(roots_c), 288);
	roots_c = roots_c * -1;
	EXPECT_EQ(double(roots_c), -288);
}
