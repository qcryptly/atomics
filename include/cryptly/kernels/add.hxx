#pragma once
#include "cryptly/cuda_ops.hxx"

// Stopper function

template<typename ...>
struct add_impl {
	int operator()(){};
};
// Define and declare all possible add_impl types
template<typename T, typename ... Ts>
struct add_impl<T, Ts...> : public add_impl<Ts...> {
	using add_impl<Ts...>::operator();
	int operator()(const kernel_ps_t&,
		std::size_t,
		const T*,
		const T*, T*);
};


struct add : public add_impl<primitive_types> {};
