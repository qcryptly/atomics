#pragma once

#include <boost/rational.hpp>
#include <ostream>

#include "cryptly/math/rational.hxx"
#include "cryptly/math/linear_ops.hxx"

namespace cryptly::math {
inline namespace v1 {

using rational_t = std::int64_t;

template<class TArg>
struct jm_s {
	using unit_t = cryptly::math::brational_t<TArg>;
	unit_t j{};
	unit_t m{};
};

template<class TArg>
bool operator==(const jm_s<TArg>& a, const jm_s<TArg>& b){
	return a.j == b.j && a.m == b.m;
}	

template<class TArg>
std::size_t coupled_hash(const std::tuple<jm_s<TArg>, jm_s<TArg>>& state) {
	auto[a, b] = state;
	return std::hash<double>{}(boost::rational_cast<double>(a.j)) ^
		   std::hash<double>{}(boost::rational_cast<double>(a.m)) ^
		   std::hash<double>{}(boost::rational_cast<double>(b.j)) ^
		   std::hash<double>{}(boost::rational_cast<double>(b.m));
}

template<class TArg>
constexpr auto srt_coupled = [](const TArg& a, const TArg& b){
	auto [a_ff, a_fs] = a.first;
	auto [b_ff, b_fs] = b.first;
	if (a_ff.m > b_ff.m) return true;
	if (a_ff.m == b_ff.m) {
		if (a_fs.m >= b_fs.m) return true;
	}
	return false;
};

template<class TArg = std::int64_t>
struct angular_s {
	using jm_t = jm_s<TArg>;
	using unit_t = cryptly::math::brational_t<TArg>;
	using eigenvalue_t = cryptly::math::rational_s<TArg>;
	using component_t = typename cryptly::math::rational_s<TArg>::component_t;
	using state_t = std::tuple<jm_t, jm_t>;
	using state_hash_t = std::unordered_map<
		state_t, eigenvalue_t, decltype(&coupled_hash<TArg>)>;
	using spectrum_t = cryptly::math::matrix_s<state_t, eigenvalue_t>;
	using par_t = typename spectrum_t::par_t; 

	spectrum_t _spectrum{};
	const jm_t _total{};

	double eps = cryptly::math::high_eps;

	angular_s(spectrum_t& spectrum, jm_t total) : 
		_spectrum{spectrum}, _total{total} {
		_spectrum.eps = eps;
	};

	void set_eps(double _eps) {
		_spectrum.eps = eps;
		eps = _eps;
	}

	friend std::ostream &operator<<(std::ostream& os, angular_s& ang) {
		jm_t current = ang._total;
		for (auto wf : ang._spectrum._matrix) {
			os << "J = " << current.j << ", M = " << current.m << "\n";
			current.j--;
			for (auto element : wf) {
				auto[jm_a, jm_b] = element.first;
				os << "(" << double(element.second) << ")"
					<< "|" << jm_a.j << "," << jm_a.m << ">"
					<< "|" << jm_b.j << "," << jm_b.m << ">"
					<< "\n";
			}
		}
		return os;
	}	
	eigenvalue_t down(const jm_t& j) const {
		return eigenvalue_t{
			component_t{
				// scalar
				{1},
				// root
				{(j.j - j.m + 1)*(j.j + j.m)},
				// power
				{1,2}
			}};
	};

	eigenvalue_t get_ev(state_hash_t& state_hash,
			const state_t& side,
			const jm_t& j,
			const eigenvalue_t& scalar,
			const eigenvalue_t& norm_const) const {
		if (state_hash.find(side) == state_hash.end()) {
			// We found ev
			return (down(j) * scalar) * norm_const;
		} else {
			return ((down(j) * scalar) * norm_const) + state_hash[side];
		}
	};

	angular_s jdown(){
		// Now that we have a wf lets build a new
		// wf
		spectrum_t spectrum{srt_coupled<par_t>};
		spectrum.eps = eps;
		auto _j = _total.j;
		auto _m = _total.m;
		for(auto wf : _spectrum._matrix) {
			const auto norm_const = down(jm_t{_j, _m}).invert();
			_j--;
			state_hash_t new_state{{}, coupled_hash<TArg>};
			for (auto state : wf) {
				auto[jm_a, jm_b] = state.first;
				const auto scalar = state.second;
				state_t left {{jm_a.j, jm_a.m - 1}, jm_b};
				state_t right {jm_a, {jm_b.j, jm_b.m - 1}};
				new_state[left] = get_ev(new_state, left, jm_a, scalar, norm_const);
				new_state[right] = get_ev(new_state, right, jm_b, scalar, norm_const);
			}
			spectrum.push_back(new_state);
		}
		spectrum.push_back(spectrum.cross());
		auto result = angular_s{spectrum, jm_t{_total.j, _total.m - 1}};
		result.set_eps(eps);
		return result;
	};
};

}
}
