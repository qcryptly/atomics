#include "cryptly/test.hxx"

#include "cryptly/math/linear_ops.hxx"
#include "cryptly/math/wavefunction.hxx"

using number_t = double;
using angular_t = typename cryptly::math::angular_s<number_t>;
using jm_t = typename cryptly::math::angular_s<number_t>::jm_t;
using eigenvalue_t = typename angular_t::eigenvalue_t;
using state_t = typename angular_t::state_t;
using matrix_state_t =
    typename cryptly::math::matrix_s<state_t, eigenvalue_t>::state_t;
using matrix_states_t =
    typename cryptly::math::matrix_s<state_t, eigenvalue_t>::matrix_t;
using state_hash_t = angular_t::state_hash_t;
using par_t = typename cryptly::math::matrix_s<state_t, eigenvalue_t>::par_t;
using matrix_t = cryptly::math::matrix_s<state_t, eigenvalue_t>;

TEST(wavefunction, jdown) {
  constexpr auto st = 5.0 / 2.0;
  jm_t maximal{5, 5};

  // Set our initial key state
  state_t initial_state{{st, st},  // |(5/2),(5/2)>
                        {st, st}}; // |(5/2),(5/2)>

  // Create initial wf vector
  matrix_state_t initial_wf{{initial_state, 1.0}};
  matrix_t spectrum{cryptly::math::srt_coupled<par_t>};
  spectrum.push_back(initial_wf);

  angular_t angular_wf{spectrum, maximal};
  angular_wf.set_eps(cryptly::math::med_eps);
  std::vector<angular_t> states{angular_wf};
  for (auto i = 5.0; i > 0.0; i--)
    states.push_back(states[states.size() - 1].jdown());

  angular_t angular_wf_expected_0{spectrum, maximal};
  EXPECT_EQ(states[0], angular_wf_expected_0);

  double half_rt = std::sqrt(1.0 / 2.0);
  angular_t angular_wf_expected_1{
      matrix_t{matrix_states_t{
                   // First wf vector
                   {{{{2.5, 2.5}, {2.5, 1.5}}, half_rt},
                    {{{2.5, 1.5}, {2.5, 2.5}}, half_rt}},
                   {{{{2.5, 2.5}, {2.5, 1.5}}, half_rt},
                    {{{2.5, 1.5}, {2.5, 2.5}}, -1.0 * half_rt}},
               },
               cryptly::math::srt_coupled<par_t>},
      maximal};
  EXPECT_EQ(states[1], angular_wf_expected_1);
  for (auto x : states)
    std::cout << x << std::endl;
}
