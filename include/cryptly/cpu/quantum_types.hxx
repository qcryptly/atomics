#include <boost/rational.hpp>

namespace cryptly {
inline namespace v1 {

using rational_t = std::int64_t;
using principal_t = boost::rational<rational_t>;
using orbital_t = boost::rational<rational_t>;
using magnetic_t = boost::rational<rational_t>;
using spin_t = boost::rational<rational_t>; 

using j_t = boost::rational<rational_t>;
using m_t = boost::rational<rational_t>;

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

}
}
