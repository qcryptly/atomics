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

ev_vect_t operator*(const ev_vect_t& vec_a, const ev_vect_t& vec_b) {
	ev_vect_t results{};
	for (auto p_a : vec_a) {
		for(auto p_b : vec_b) {
			const auto key = p_a.first * p_b.first;
			const auto scalar = p_a.second * p_b.second;
			if (results.find(key) == results.end()) {
				results[key] = scalar;
				continue;
			}
			results[key] = results[key] * scalar;
		}
	}
	return results;
}

ev_vect_t operator+(const ev_vect_t& vec_a, const ev_vect_t& vec_b) {
	ev_vect_t results{};
	// Add the first vector
	for (auto e : vec_a) {
		auto key = e.first;
		auto scalar = e.second;
		if (results.find(key) == results.end()) {
			results[key] = scalar;
			continue;
		}
		results[key] += scalar;
	}
	// Add the second vector
	for (auto e : vec_b) {
		auto key = e.first;
		auto scalar = e.second;
		if (results.find(key) == results.end()) {
			results[key] = scalar;
			continue;
		}
		results[key] += scalar;
	}
	return results;
}

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

constexpr auto print_spectrum = [](const auto& spectrum, auto& log){
	for (auto wf : spectrum) print_wf(wf, log);
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
constexpr auto j_down = [](const wfj_t& _wfj, const j_t J, const m_t M){
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

template<class TKey, class TVal>
struct matrix_s {
	using key_t = TKey;
	using val_t = TVal;

	using state_t = std::vector<std::pair<TKey, TVal>>;
	std::vector<state_t> _matrix{};
	matrix_s() {};
	~matrix_s(){};
	void push_back(std::vector<std::pair<TKey, TVal>> mvalues){
		_matrix.push_back(mvalues);
	}
	void push_back(std::pair<TKey, TVal> mvalue){
		_matrix.push_back(mvalue);
	}
	const auto size() {
		return _matrix.size();
	}
	// Not thread safe
	state_t &operator[](std::size_t i){
		return _matrix[i];
	}
};

// Expect matrix = (n)(n)
template <class TArg = std::int64_t>
constexpr auto det = [](auto& matrix){
	auto size = matrix.size();
	TArg agg{};
	for(auto col = 0; col < size; col ++) {
		TArg element = col % 2 == 0 ? 1 : -1;
		for(auto row = 0; row < size; row ++) {
			// matrix[][].second == ev_vect_t
			// we need to broadcast scalars to the element
			//
			// we also need the element to be of type ev_vect_t
			element = (element * matrix[row][(row + col) % size].second);
		}
		agg = agg + element;
	}
	return agg;	
};

// Expect matrix = (n - 1) = (m)
template<class TArg = std::int64_t>
constexpr auto cross = [](auto& matrix, auto& orth_vect){
	auto size = matrix.size();
	TArg agg{};
	using __matrix_t = typename std::remove_reference<decltype(matrix)>::type;
	using __vector_t = typename std::remove_reference<decltype(orth_vect)>::type;
	for(auto col = 0; col < size + 1; col ++) {
		auto element = col % 2 == 0 ? 1 : -1;
		__matrix_t sub_matrix{};
		// Create sub matrix
		for (auto row = 0; row < size; row ++) {
			__vector_t new_vect{};
			for(auto col_i = 0; col_i < size + 1; col_i ++) {
				if (col_i != col)
					new_vect.push_back(matrix[row][col_i]);
			}
			// Uses copy constructor, this is slow
			sub_matrix.push_back(new_vect);
		}
		// Set our vector
		orth_vect[col].second = det<TArg>(sub_matrix) * element;
	}
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
