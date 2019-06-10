#include "cryptly/test.hxx"

#include "cryptly/math/rational.hxx"
#include "cryptly/math/linear_ops.hxx"
#include "cryptly/math/wavefunction.hxx"

using number_t = std::int64_t;
using rational_t = typename cryptly::math::rational_s<number_t>;
using angular_t = typename cryptly::math::angular_s<number_t>;
using jm_t = typename cryptly::math::angular_s<number_t>::jm_t;
using eigenvalue_t = typename angular_t::eigenvalue_t;
using state_t = typename angular_t::state_t;
using hashmap_t = typename cryptly::math::rational_s<number_t>::hashmap_t;
using hashkey_t = typename cryptly::math::rational_s<number_t>::hashkey_t;
using matrix_state_t = typename cryptly::math::matrix_s<state_t, eigenvalue_t>::state_t;
using state_hash_t = angular_t::state_hash_t;
using par_t = typename cryptly::math::matrix_s<state_t, eigenvalue_t>::par_t;
using matrix_t = cryptly::math::matrix_s<state_t, eigenvalue_t>;

TEST(wavefunction, jdown) {
	jm_t maximal{5,5};

	// Set our initial key state
	state_t initial_state{
		{{5,2},{5,2}}, // |(5/2),(5/2)>
		{{5,2},{5,2}}};// |(5/2),(5/2)>

	// Set our initial eigenvalue
	hashmap_t one_hash{{
		{{1,1},1}
	}, 1, cryptly::math::hashkey_rational_h<hashkey_t>};	
	eigenvalue_t initial_eigenvalue{one_hash};

	// Create initial wf vector	
	matrix_state_t initial_wf{{initial_state, initial_eigenvalue}};	
	matrix_t spectrum{cryptly::math::srt_coupled<par_t>};
	spectrum.push_back(initial_wf);

	angular_t angular_wf{spectrum, maximal};
	angular_wf.set_eps(cryptly::math::med_eps); 
	std::vector<angular_t> states{angular_wf};
	for(auto i = 5; i > 0; i--) states.push_back(states[states.size() - 1].jdown());
	for(auto x : states) std::cout << x << std::endl;

}
