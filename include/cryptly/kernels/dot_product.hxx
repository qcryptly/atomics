#pragma once
#include "cryptly/cuda_ops.hxx"

namespace cryptly {
inline namespace v1 {

template<typename...>
struct dot_product_impl {
	void operator()(){};
};

template<typename T, typename...Ts>
struct dot_product_impl<T, Ts...> : public dot_product_impl<Ts...> {
	using dot_product_impl<Ts...>::operator();

	int operator()(const kernel_ps_t&,
		std::size_t, const T*, const T*, T*);
};

struct dot_product : public dot_product_impl<primitive_types>{};

}
}
