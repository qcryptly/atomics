#pragma once

#include <boost/rational.hpp>
#include <cmath>
#include <ostream>
#include <limits>
#include <type_traits>

namespace cryptly::math {
inline namespace v1 {

enum class component_e : std::int32_t {scalar, root, power};

template<class TArg = std::int64_t>
using brational_t = boost::rational<TArg>;

// scalar, root, power
template<class TArg = std::int64_t>
using rcomponent_t = std::tuple<
	brational_t<TArg>, brational_t<TArg>, brational_t<TArg>>;

constexpr auto scalar = [](auto input){
	return std::get<static_cast<std::int32_t>(component_e::scalar)>(input);
};

constexpr auto root = [](auto input){
	return std::get<static_cast<std::int32_t>(component_e::root)>(input);
};

constexpr auto power = [](auto input){
	return std::get<static_cast<std::int32_t>(component_e::power)>(input);
};

constexpr auto evaluate = [](auto component){
	return std::pow(
		boost::rational_cast<double>(root(component)),
		boost::rational_cast<double>(power(component))) *
		boost::rational_cast<double>(scalar(component));
};

template<class TArg, typename>
struct rational_s;

template<class TArg>
std::size_t hashkey_rational_h(const TArg& key) {
	return std::hash<double>{}(std::pow(
			boost::rational_cast<double>(key.first),
			boost::rational_cast<double>(key.second)));
}

template<class TComponent>
constexpr auto to_component = [](const auto& input){
	return TComponent{input.second,
					  input.first.first,
					  input.first.second};
};

template<class TComponent>
constexpr auto set_hash_to_vect = [](const auto& hm, auto& vector_num){
	vector_num.clear();
	for (auto element : hm) {
		vector_num.push_back(to_component<TComponent>(element));
	}
	std::sort(vector_num.begin(), vector_num.end(), [](auto comp_a, auto comp_b) {
		return evaluate(comp_a) > evaluate(comp_b);	
	});
};

template<class TKey>
constexpr auto get_key = [](const auto& comp){
	return TKey{root(comp), power(comp)};
};

template<class TKey, class TVal>
constexpr auto set_vect_to_hash = [](const auto& vect, auto& hashmap){
	hashmap.clear();
	for (auto element : vect) {
		auto _scalar = scalar(element);
		auto key = get_key<TKey>(element);
		if (hashmap.find(key) == hashmap.end()) {
			hashmap[key] = _scalar;
			continue;
		}
		hashmap[key] += _scalar;
	}
};

template<class TKey>
constexpr auto add_hash = [](auto& acc, const auto& input) {
	auto key = get_key<TKey>(input);
	// copy our hash num
	if (acc.find(key) == acc.end()) {
		acc[key] = scalar(input);
	} else {
		acc[key] += scalar(input);
	}
};

template<class TSize>
constexpr auto to_rational = [](const double input) {
	TSize num{};
	TSize dem{};
	if (std::abs(input) <= 1) {
		dem = std::numeric_limits<TSize>::max();
		num = TSize(double(dem) * input);
	} else {
		num = std::numeric_limits<TSize>::max();
		dem = TSize(double(1.0/input) * double(num));
	}
	return boost::rational<TSize>{num, dem};
};
// Unfortunately, I don't know how to generically
// convert our expression into a single rational radical
// w/ generic power. I'm sure there's something in abstract
// algebra that can tell me if it's possible or not,
// but multiplication will have to kill percision at this point
template<class TSize>
constexpr auto mul_vect = [](auto acc, const auto& input) {
	for (auto& element : acc) {
		std::get<0>(element) = to_rational<TSize>(evaluate(input) * evaluate(element));
		std::get<1>(element) = 1;
		std::get<2>(element) = 1;
	}
	return acc;
};

template<class TArg, typename = std::enable_if_t<std::is_integral_v<TArg>, TArg>>
struct rational_s {
	using primitive_t = TArg;
	using base_t = brational_t<TArg>;
	using component_t = rcomponent_t<TArg>;
	using vector_t = std::vector<rcomponent_t<TArg>>;
	// Use the root and power to form a key
	using hashkey_t = std::pair<base_t, base_t>;
	using hashval_t = base_t;
	using hash_t = decltype(&hashkey_rational_h<hashkey_t>);
	// key = <root, power>, value = scalar
	using hashmap_t = std::unordered_map<
		hashkey_t, hashval_t, hash_t>;

	vector_t vector_num{};
	hashmap_t hash_num{{}, hashkey_rational_h<hashkey_t>};

	rational_s(hashmap_t hm) : hash_num{hm} {
		set_hash_to_vect<component_t>(hm, vector_num);
	};
 
	rational_s(vector_t vect) : vector_num{vect} {
		set_vect_to_hash<hashkey_t, hashval_t>(vect, hash_num);
	};

	rational_s(component_t element) {
		vector_num.push_back(element);	
		set_vect_to_hash<hashkey_t, hashval_t>(vector_num, hash_num);
	};

	template<typename TA, std::enable_if_t<std::is_same_v<TA, base_t>, base_t>>
	explicit rational_s(TA element) {
		vector_num.push_back(component_t{
			element, // scalar
			base_t{1,1},	
			base_t{1,1}});	
		set_vect_to_hash<hashkey_t, hashval_t>(vector_num, hash_num);
	};

	template<typename TA>
	rational_s(TA element) {
		if constexpr(std::is_same_v<TA, double>) {
			vector_num.push_back(component_t{
				to_rational<TArg>(element), // scalar
				base_t{1,1},	
				base_t{1,1}});	
			set_vect_to_hash<hashkey_t, hashval_t>(vector_num, hash_num);
		} else if(std::is_same_v<TA, primitive_t>) {
			vector_num.push_back(component_t{
				base_t{element,1}, // scalar
				base_t{1,1},	
				base_t{1,1}});	
			set_vect_to_hash<hashkey_t, hashval_t>(vector_num, hash_num);
		}
	};

	rational_s() : vector_num{}, hash_num{} {};

	explicit operator double() const noexcept{
		double agg{};
		for (auto element : vector_num) {
			agg += evaluate(element);
		}
		return agg;
	}

	friend std::ostream &operator<<(std::ostream& os, const rational_s& rn){
		auto count = 0;
		for (auto element : rn.vector_num) {
			auto _scalar = scalar(element);
			auto _root = root(element);
			auto _power = power(element);
			os << "(" << _scalar << ")root(" << _root
			<< ")^(" << _power << ")";
			if (count < rn.vector_num.size() -1) {
				os << " + ";
			}
			count++;
		}
		return os;
	}	

	bool operator==(const rational_s& rn) const noexcept {
		return std::abs(double((*this)) -
					double(rn)) < std::numeric_limits<double>::epsilon();
	}

	bool operator!=(const rational_s& rn) const noexcept {
		return not (*this == rn);
	}

	bool operator>(const rational_s& rn) const noexcept {
		return double((*this)) > double(rn);
	}

	rational_s invert() {
		vector_t new_vect{};
		for(auto element : vector_num) {
			auto _scalar = scalar(element);
			auto _root = root(element);
			auto _power = power(element);
			new_vect.push_back(component_t{
				1 / _scalar,
				_root,
				-1 * _power
			});
		}
		return rational_s<TArg>{new_vect};
	}

	/* Addition */
	rational_s operator+(const rational_s& input) {
		return *this + input.hash_num;	
	}
	rational_s operator+(const vector_t& input) {
		hashmap_t new_hash_num{};
		set_vect_to_hash<hashkey_t, hashval_t>(input, new_hash_num);
		return *this + new_hash_num;	
	}
	rational_s operator+(const hashmap_t& input) {
		hashmap_t new_hash{hash_num.begin(), hash_num.end(), input.size(), hashkey_rational_h<hashkey_t>};
		for(const auto element : input) {
			add_hash<hashkey_t>(new_hash, to_component<component_t>(element));
		}
		return rational_s<TArg>{new_hash};	
	}
	rational_s operator+(const component_t& input) {
		hashmap_t new_hash{hash_num.begin(), hash_num.end()};
		add_hash<hashkey_t>(new_hash, input);
		return rational_s<TArg>{new_hash};
	}
	rational_s operator+(const base_t& input) {
		rational_s<TArg> new_number{};
		
		return *this + component_t{
			input,
			base_t{1,1},
			base_t{1,1}	
		};	
	}
	rational_s operator+(const TArg& input) {
		rational_s<TArg> new_number{};
		return *this + component_t{
			base_t{input,1}, // scalar
			base_t{1,1},	
			base_t{1,1}};	
	}

	/* Multiplication */
	rational_s operator*(const rational_s& input) {
		return *this * input.vector_num;	
	}
	rational_s operator*(const vector_t& input) {
		vector_t new_vect{
			vector_num.begin(), vector_num.end()};
		vector_t result{};
		// This effectivel converts all roots and powers to (1)(1)
		// with varying scalars...
		for(const auto element : input) {
			auto r = mul_vect<TArg>(new_vect, element);
			result.insert(result.end(), r.begin(), r.end());
		}
		return rational_s<TArg>{result};	
	}
	rational_s operator*(const hashmap_t& input) {
		vector_t new_vect{};
		set_hash_to_vect<component_t>(input, new_vect);
		return *this * new_vect;	
	}
	rational_s operator*(const component_t& input) {
		return rational_s<TArg>{mul_vect<TArg>(vector_num, input)};
	}
	rational_s operator*(const base_t& input) {
		return *this * component_t{
			input,
			base_t{1,1},
			base_t{1,1}	
		};	
	}
	rational_s operator*(const primitive_t& input) {
		return *this * component_t{
			base_t{input,1}, // scalar
			base_t{1,1},	
			base_t{1,1}};	
	}
}; 

}
}
