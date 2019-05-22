#pragma once

#include <boost/rational.hpp>
#include <tuple>
#include <vector>
#include <array>
#include <unordered_map>

namespace cryptly {
inline namespace v1 {

using rational_t = std::int64_t;
using principal_t = boost::rational<rational_t>;
using orbital_t = boost::rational<rational_t>;
using magnetic_t = boost::rational<rational_t>;
using spin_t = boost::rational<rational_t>; 

using j_t = boost::rational<rational_t>;
using m_t = boost::rational<rational_t>;
using atomic_t = boost::rational<rational_t>;
using ev_t = boost::rational<rational_t>;


struct state_s {
	principal_t n{};	
	orbital_t l{};	
	magnetic_t m{};	
	spin_t x{};	
};

struct jm_s {
	j_t j{};
	m_t m{};
};

bool operator==(const jm_s& a, const jm_s& b){
	return a.j == b.j && a.m == b.m;
}	

using wfj_tup_t = std::tuple<jm_s, jm_s>;

std::size_t wfj_hash(const wfj_tup_t& wfj) {
	wfj_tup_t x = wfj;
	auto[a, b] = x;
	return std::hash<double>{}(boost::rational_cast<double>(a.j)) ^
		   std::hash<double>{}(boost::rational_cast<double>(a.m)) ^
		   std::hash<double>{}(boost::rational_cast<double>(b.j)) ^
		   std::hash<double>{}(boost::rational_cast<double>(b.m));
}

std::size_t rational_hash(const ev_t& input) {
	return std::hash<double>{}(
		double(input.numerator()) /
		double(input.denominator())
	);
}

using ev_vect_t = std::unordered_map<ev_t, rational_t, decltype(&rational_hash)>;
using wfj_t = std::unordered_map<wfj_tup_t, ev_vect_t, decltype(&wfj_hash)>;

constexpr auto print_ev_vect = [](const auto& ev_vect){
	double result = {};
	for (auto& ev : ev_vect) {
		result += double(ev.second) * std::sqrt(boost::rational_cast<double>(ev.first));
	}
	return result;
};

constexpr auto print_wf = [](const auto& wf, auto& log){
	for (auto element : wf) {
		auto[jm_a, jm_b] = element.first;
		log << "(" << print_ev_vect(element.second) << ")"
			<< "|" << jm_a.j << "," << jm_a.m << ">"
			<< "|" << jm_b.j << "," << jm_b.m << ">"
			<< "\n";
	}
};

constexpr auto get_ev_down = [](const jm_s& j){
	return (j.j - j.m + 1)*(j.j + j.m);
};

constexpr auto get_total_norm_const = [](const j_t& j, const m_t& m){
	jm_s init{j,m};
	return get_ev_down(init);
};

constexpr auto get_root_hash = [](
	const auto& jm,
	const auto& root_key,
	const auto& norm_const){
	return (get_ev_down(jm) * root_key) / norm_const;
};

constexpr auto set_ev = [](auto& wfj,
		const auto& root_key, const auto& scalar,
		const auto& side, const auto& jm,
		const auto& norm_const){
	const auto root_hash = get_root_hash(jm, root_key, norm_const);
	if (wfj.find(side) == wfj.end()) {
		// We found ev
		ev_vect_t init_ev {{}, cryptly::rational_hash};
		init_ev[root_hash] = scalar;
		wfj[side] = init_ev;
	} else {
		// Will give us a hash of our root
		if (wfj[side].find(root_hash) == wfj[side].end()) {
			wfj[side][root_hash] = scalar;
		} else {
			wfj[side][root_hash] += scalar;
		}
	}
};
constexpr auto j_down = [](const wfj_t& _wfj, j_t J, m_t M){
	// Now that we have a wf lets build a new
	// wf
	const auto norm_const = get_total_norm_const(J, M);
	wfj_t wfj{{}, wfj_hash};
	for (auto tup_e : _wfj) {
		auto[jm_a, jm_b] = tup_e.first;
		for (auto &ev : tup_e.second) {
			auto root_hash = ev.first;
			auto scalar = ev.second;
			// For each state perform a ladder
			wfj_tup_t left {{jm_a.j, jm_a.m - 1}, jm_b};
			wfj_tup_t right {jm_a, {jm_b.j, jm_b.m - 1}};
			set_ev(wfj, root_hash, scalar, left, jm_a, norm_const);
			set_ev(wfj, root_hash, scalar, right, jm_b, norm_const);
		}
	}
	
	return wfj;
};

template<class TArg = double>
struct cg final {
	jm_s j_a{};
	jm_s j_b{};
	TArg operator()(jm_s total) noexcept {
		return std::sqrt(1.0 / 2.0);
	};
};


}
}
