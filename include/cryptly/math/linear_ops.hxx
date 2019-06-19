#pragma once

#include <cassert>
#include <cmath>
#include <type_traits>
#include <vector>

namespace cryptly::math {
inline namespace v1 {

template <class TArg>
constexpr auto ascending =
    [](const TArg &a, const TArg &b) { return !(a.first > b.first); };

constexpr double coarse_eps = 1.0e-3;
constexpr double bad_eps = 1.0e-6;
constexpr double low_eps = 1.0e-9;
constexpr double med_eps = 1.0e-14;
constexpr double high_eps = std::numeric_limits<double>::epsilon();

template <class TKey, class TVal> struct matrix_s {
  using key_t = TKey;
  using val_t = TVal;
  using par_t = std::pair<TKey, TVal>;
  using state_t = std::vector<std::pair<TKey, TVal>>;
  using matrix_t = std::vector<state_t>;
  matrix_t _matrix{};
  bool (*_srt)(const par_t &, const par_t &);
  TVal eps = high_eps;
  matrix_s(bool (*srt)(const par_t &, const par_t &)) { _srt = srt; };
  matrix_s(const matrix_t matrix) : _matrix{matrix} {
    _srt = ascending<par_t>;
  };
  matrix_s(const matrix_t matrix, bool (*srt)(const par_t &, const par_t &))
      : _matrix{matrix} {
    _srt = srt;
  };
  matrix_s() { _srt = ascending<par_t>; };
  ~matrix_s(){};
  void push_back(state_t mvalues) { _matrix.push_back(mvalues); }
  template <class... TArg>
  void push_back(std::unordered_map<TArg...> hash_num) {
    // Must convert to vector
    state_t new_vector{};
    for (auto element : hash_num) {
      new_vector.push_back(par_t{element.first, element.second});
    }
    _matrix.push_back(new_vector);
  }
  void push_back(std::pair<TKey, TVal> mvalue) { _matrix.push_back(mvalue); }
  const auto size() { return _matrix.size(); }
  // Not thread safe
  state_t &operator[](std::size_t i) { return _matrix[i]; }

  bool operator==(const matrix_s &matrix) const {
    return *this == matrix._matrix;
  }
  bool operator==(const matrix_t &matrix) const {
    if (_matrix.size() != matrix.size())
      return false;
    for (auto row_index = 0; row_index < _matrix.size(); ++row_index) {
      auto row_left = _matrix[row_index];
      auto row_right = matrix[row_index];
      if (row_left.size() != row_right.size())
        return false;
      for (auto col_index = 0; col_index < row_left.size(); ++col_index) {
        auto val_left = row_left[col_index].second;
        auto val_right = row_right[col_index].second;
        auto key_left = row_left[col_index].first;
        auto key_right = row_right[col_index].first;
        auto diff = std::abs(val_left - val_right);
        if (diff > eps or key_left != key_right) {
          return false;
        }
      }
    }
    return true;
  }

  TVal det(matrix_s<TKey, TVal> &matrix) {
    auto size = matrix.size();
    TVal agg{};
    for (auto col = 0; col < size; col++) {
      TVal element = TVal(col % 2 == 0 ? 1 : -1);
      for (auto row = 0; row < size; row++) {
        // matrix[][].second == ev_vect_t
        // we need to broadcast scalars to the element
        //
        // we also need the element to be of type ev_vect_t
        element = (element * matrix[row][(row + col) % size].second);
      }
      agg = agg + element;
    }
    return agg;
  }
  // orth_vector must be initialized with the correct ordering of states
  void cross(state_t &orth_vect) {
    auto size = _matrix.size();
    TVal agg{};
    for (auto col = 0; col < size + 1; col++) {
      TVal element = TVal(col % 2 == 0 ? 1 : -1);
      matrix_s<TKey, TVal> sub_matrix{_srt};
      // Create sub matrix
      for (auto row = 0; row < size; row++) {
        state_t new_vect{};
        for (auto col_i = 0; col_i < size + 1; col_i++) {
          if (col_i != col)
            new_vect.push_back(_matrix[row][col_i]);
        }
        // Uses copy constructor, this is slow
        sub_matrix.push_back(new_vect);
      }
      // Set our vector
      orth_vect[col].second = det(sub_matrix) * element;
    }
  }

  state_t mul(state_t vect, TVal value) {
    for (auto &el : vect) {
      el.second *= value;
      if (std::abs(el.second) < eps) {
        el.second = 0;
      }
    }
    return vect;
  }

  state_t add(state_t vect_a, const state_t &vect_b) {
    for (int i = 0; i < vect_a.size(); ++i) {
      auto &num = vect_a[i].second;
      num += vect_b[i].second;
      if (std::abs(num) < eps) {
        num = 0;
      }
    }
    return vect_a;
  }

  // Creates copy of matrix... still immutable :)
  matrix_t reduce_upper(matrix_t matrix) {
    // Row echelon
    auto row_size = 0;
    if (matrix.size() > 0)
      row_size = matrix[0].size();
    for (auto row = matrix.end() - 1; row != matrix.begin() - 1; --row) {
      auto index = 0;
      // Zero out leading coes
      auto coe = (*row)[index].second;
      while (std::abs(coe - 1.0) > eps) {
        index++;
        coe = (*row)[index].second;
        if (index + 1 > row_size)
          continue;
      }

      for (auto row_i = row - 1; row_i != matrix.begin() - 1; --row_i) {
        auto coe_i = (*row_i)[index].second;
        auto r = mul(*row, coe_i * -1.0);
        *row_i = add(*row_i, r);
      }
    }
    return matrix;
  }

  // Creates copy of matrix... still immutable :)
  matrix_t zero_upper(matrix_t matrix) {
    auto index = 0;
    auto row_size = 0;
    if (matrix.size() > 0)
      row_size = matrix[0].size();
    // Row echelon
    for (auto row = matrix.begin(); row != matrix.end() and index < row_size;
         ++row) {
      auto coe = (*row)[index].second;
      auto tmp = *row;
      // Sort leading zero
      auto row_i = row + 1;
      while (std::abs(coe) < eps) {
        if (matrix.end() == row_i) {
          row_i = row + 1;
          if (index + 1 < row_size) {
            index++;
            coe = (*row)[index].second;
            continue;
          } else {
            return matrix;
          }
        } else {
          auto coe_i = (*row_i)[index].second;
          auto current = *row_i;
          if (std::abs(coe_i) > eps) {
            *row = *row_i;
            *row_i = tmp;
            coe = coe_i;
            break;
          }
          row_i++;
        }
      }

      *row = mul(*row, 1.0 / coe);
      // Zero out leading coes
      for (row_i = row + 1; row_i != matrix.end(); ++row_i) {
        auto coe_i = (*row_i)[index].second;
        auto r = mul(*row, -1.0 * coe_i);
        *row_i = add(*row_i, r);
      }
      ++index;
    }
    return matrix;
  }

  matrix_t rref() {
    // Copy
    matrix_t matrix = zero_upper(_matrix);
    matrix = reduce_upper(matrix);
    // Reduce leads to zero
    return matrix;
  }

  TVal dot(state_t vec_a, state_t vec_b) {
    assert(vec_b.size() == vec_a.size());
    TVal agg{};
    for (auto index = 0; index < vec_b.size(); index++) {
      agg += vec_a[index].second * vec_b[index].second;
    }
    return agg;
  }

  state_t cross() {
    state_t vect_orth{};
    std::vector<bool> is_set{};

    auto matrix = rref();

    auto row_size = 0;
    if (matrix.size() > 0) {
      auto row = matrix[0];
      row_size = row.size();
      is_set.resize(row_size);
      vect_orth.resize(row_size);
      std::fill(is_set.begin(), is_set.end(), false);
      // We don't care what values are, they will be overwritten
      // we only care about preserving key states
      std::copy(row.begin(), row.end(), vect_orth.begin());
    }

    for (auto row = matrix.end() - 1; row != matrix.begin() - 1; row--) {
      auto index = 0;
      auto leading_index = -1;
      // on the first row
      for (auto element : *row) {
        if (leading_index > -1) {
          if (not is_set[index]) {
            vect_orth[index].second = 1;
            vect_orth[leading_index].second -= element.second;
            is_set[index] = true;
          } else {
            vect_orth[leading_index].second -=
                vect_orth[index].second * element.second;
          }
        } else if (std::abs(element.second - 1.0) < eps) {
          leading_index = index;
          vect_orth[index].second = 0.0;
          is_set[leading_index] = true;
        }
        ++index;
      }
    }

    auto metric = 1.0 * dot(vect_orth, vect_orth);

    metric = (matrix.size() % 2 == 0 ? 1.0 : -1.0) / std::sqrt(metric);
    vect_orth = mul(vect_orth, metric);

    return vect_orth;
  }
};

} // namespace v1
} // namespace cryptly::math
