#pragma once
#include <iostream>

#define primitive_types std::uint8_t, std::uint16_t, \
						std::uint32_t, std::uint64_t, \
						std::int8_t, std::int16_t, \
						std::int32_t, std::int64_t, \
						float, double, bool

auto peak_mem_band = [](auto bus_width, auto clock_rate) -> float {
	return 2.0*clock_rate*(bus_width/8)/1.0e+6;
};

constexpr auto cuda_success{0};


int print_devices();

constexpr auto default_grid_size = 1;
constexpr auto default_block_size = 1;
constexpr auto default_shared_size = 48 * (1 << 10);
constexpr auto default_cuda_stream = 0;
struct kernel_ps_t {
	std::size_t grid_size{default_grid_size};
	std::size_t block_size{default_block_size};
	std::size_t shared_size{default_shared_size};
	std::size_t cuda_stream{default_cuda_stream};
};
