//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_ARRAY_HH_
#define ELM_KERNEL_TEST_ARRAY_HH_

#include <algorithm>
#include <array>
#include <assert.h>
#include <memory>
#include <numeric>
#include <string>

namespace ELM {
namespace Impl {

template <typename T> class Data_ {
public:
  // assigment operator
  void operator=(T t) {
    for (auto &d : *this) {
      d = t;
    }
  }

  // iterators
  T const *begin() const { return &(d_[0]); }
  T const *end() const { return &(d_[len_]); }
  T *begin() { return &(d_[0]); }
  T *end() { return &(d_[len_]); }
  size_t size() const { return len_; }
  T const *data() const { return d_; }
  T *data() { return d_; }

protected:
  //
  // Constructors are all protected -- do not use this class directly!
  //
  // construct from a total length
  Data_(int len) : len_(len), do_(std::shared_ptr<T>(new T[len], std::default_delete<T[]>())) { d_ = do_.get(); }

  // construct and initialize
  Data_(int len, T d) : Data_(len) { *this = d; }

  // construct non-owning view
  Data_(int len, T *d) : len_(len), d_(d) {}

  // copy construct
  Data_(const Data_ &other) = default;

  // resize data ptr
  // destructive
  void Data_resize(int N) {
    len_ = N;
    do_.reset(new T[N](), std::default_delete<T[]>());
    d_ = do_.get();
  }

protected:
  // global length
  int len_;

  // owning data
  std::shared_ptr<T> do_;

  // non-owning data
  T *d_;
};

} // namespace Impl

// templated class for multi-dimensional array
template <typename T, size_t D> class Array : public Impl::Data_<T> {};

// 1D specialization
template <typename T> class Array<T, 1> : public Impl::Data_<T> {
  static const size_t dim = 1;

public:
  // forward construction
  Array(int N) : Impl::Data_<T>(N) {}
  Array(const std::string& name, int N) : name_(name), Impl::Data_<T>(N) {}
  Array(std::array<int, 1> N) : Impl::Data_<T>(std::get<0>(N)) {}
  Array(const std::string& name, std::array<int, 1> N) : name_(name), Impl::Data_<T>(std::get<0>(N)) {}

  // forward construction
  Array(int N, T t) : Impl::Data_<T>(N, t) {}
  Array(const std::string& name, int N, T t) : name_(name), Impl::Data_<T>(N, t) {}
  Array(std::array<int, 1> N, T t) : Impl::Data_<T>(std::get<0>(N), t) {}
  Array(const std::string& name, std::array<int, 1> N, T t) : name_(name), Impl::Data_<T>(std::get<0>(N), t) {}

  // forward construction
  Array(int N, T *d) : Impl::Data_<T>(N, d) {}
  Array(const std::string& name, int N, T *d) : name_(name), Impl::Data_<T>(N, d) {}
  Array(std::array<int, 1> N, T *d) : Impl::Data_<T>(std::get<0>(N), d) {}
  Array(const std::string& name, std::array<int, 1> N, T *d) : name_(name), Impl::Data_<T>(std::get<0>(N), d) {}

  // forward construction
  Array(const Array<T, 1> &other) = default;

  // accessors
  T &operator()(const int i) {
    assert(0 <= i && i < len_);
    return d_[i];
  }
  const T &operator()(const int i) const {
    assert(0 <= i && i < len_);
    return d_[i];
  }

  T &operator[](const int i) {
    assert(0 <= i && i < len_);
    return d_[i];
  }
  const T &operator[](int i) const {
    assert(0 <= i && i < len_);
    return d_[i];
  }

  // resize accessor - destructive
  Array<T, 1>& resize(const int N) {
    Impl::Data_<T>::Data_resize(N);
    return *this;
  }

  // shape
  int extent(int d) const {
    assert(d < 1 && "Array::extent requested for dimension greater than is stored.");
    return shape()[d];
  }

  std::array<int, 1> shape() const { return {len_}; }

  // get variable name
  std::string label() const { return {name_}; }

  typedef T value_type;

protected:
  using Impl::Data_<T>::len_;
  using Impl::Data_<T>::d_;
  const std::string name_ = ("");
};

// 2D specialization
template <typename T> class Array<T, 2> : public Impl::Data_<T> {
  static const size_t dim = 2;

public:
  // forward construction
  Array(int M, int N) : Impl::Data_<T>(N * M), M_(M), N_(N) {}
  Array(const std::string& name, int M, int N) : name_(name), Impl::Data_<T>(N * M), M_(M), N_(N) {}
  Array(std::array<int, 2> N)
      : Impl::Data_<T>(std::get<0>(N) * std::get<1>(N)), M_(std::get<0>(N)), N_(std::get<1>(N)) {}
  Array(const std::string& name, std::array<int, 2> N)
      : name_(name), Impl::Data_<T>(std::get<0>(N) * std::get<1>(N)), M_(std::get<0>(N)), N_(std::get<1>(N)) {}

  // forward construction
  Array(int M, int N, T t) : Impl::Data_<T>(N * M, t), M_(M), N_(N) {}
  Array(const std::string& name, int M, int N, T t) : name_(name), Impl::Data_<T>(N * M, t), M_(M), N_(N) {}
  Array(std::array<int, 2> N, T t)
      : Impl::Data_<T>(std::get<0>(N) * std::get<1>(N), t), M_(std::get<0>(N)), N_(std::get<1>(N)) {}
  Array(const std::string& name, std::array<int, 2> N, T t)
      : name_(name), Impl::Data_<T>(std::get<0>(N) * std::get<1>(N), t), M_(std::get<0>(N)), N_(std::get<1>(N)) {}

  // forward construction
  Array(int M, int N, T *d) : Impl::Data_<T>(N * M, d), M_(M), N_(N) {}
  Array(const std::string& name, int M, int N, T *d) : name_(name), Impl::Data_<T>(N * M, d), M_(M), N_(N) {}
  Array(std::array<int, 2> N, T *d)
      : Impl::Data_<T>(std::get<0>(N) * std::get<1>(N), d), M_(std::get<0>(N)), N_(std::get<1>(N)) {}
  Array(const std::string& name, std::array<int, 2> N, T *d)
      : name_(name), Impl::Data_<T>(std::get<0>(N) * std::get<1>(N), d), M_(std::get<0>(N)), N_(std::get<1>(N)) {}

  // forward construction
  Array(const Array<T, 2> &other) = default;

  T &operator()(int i, int j) {
    assert(0 <= i && i < M_ && 0 <= j && j < N_);
    return d_[j + i * N_];
  }
  const T &operator()(int i, int j) const {
    assert(0 <= i && i < M_ && 0 <= j && j < N_);
    return d_[j + i * N_];
  }

  Array<T, 1> operator[](int i) {
    assert(0 <= i && i < M_);
    return Array<T, 1>(name_, N_, &d_[i * N_]);
  }
  const Array<T, 1> operator[](int i) const {
    assert(0 <= i && i < M_);
    return Array<T, 1>(name_, N_, &d_[i * N_]);
  }

  // resize accessor - destructive
  Array<T, 2>& resize(const int M, const int N) {
    Impl::Data_<T>::Data_resize(M*N);
    M_ = M;
    N_ = N;
    return *this;
  }

  // shape
  int extent(size_t d) const {
    assert(d < 2 && "Array::extent requested for dimension greater than is stored.");
    return shape()[d];
  }

  std::array<int, 2> shape() const { return {M_, N_}; }

  // get variable name
  std::string label() const { return {name_}; }

  typedef T value_type;

protected:
  int M_, N_;
  using Impl::Data_<T>::d_;
  const std::string name_ = ("");
};

// 3D specialization
template <typename T> class Array<T, 3> : public Impl::Data_<T> {
  static const size_t dim = 3;

public:
  // forward construction
  Array(int M, int N, int P) : Impl::Data_<T>(N * M * P), M_(M), N_(N), P_(P) {}
  Array(const std::string& name, int M, int N, int P) : name_(name), Impl::Data_<T>(N * M * P), M_(M), N_(N), P_(P) {}
  Array(std::array<int, 3> N)
      : Impl::Data_<T>(std::get<0>(N) * std::get<1>(N) * std::get<2>(N)), M_(std::get<0>(N)), N_(std::get<1>(N)),
        P_(std::get<2>(N)) {}
  Array(const std::string& name, std::array<int, 3> N)
      : name_(name), Impl::Data_<T>(std::get<0>(N) * std::get<1>(N) * std::get<2>(N)), M_(std::get<0>(N)), N_(std::get<1>(N)),
        P_(std::get<2>(N)) {}

  // forward construction
  Array(int M, int N, int P, T t) : Impl::Data_<T>(N * M * P, t), M_(M), N_(N), P_(P) {}
  Array(const std::string& name, int M, int N, int P, T t) : name_(name), Impl::Data_<T>(N * M * P, t), M_(M), N_(N), P_(P) {}
  Array(std::array<int, 3> N, T t)
      : Impl::Data_<T>(std::get<0>(N) * std::get<1>(N) * std::get<2>(N), t), M_(std::get<0>(N)), N_(std::get<1>(N)),
        P_(std::get<2>(N)) {}
  Array(const std::string& name, std::array<int, 3> N, T t)
      : name_(name), Impl::Data_<T>(std::get<0>(N) * std::get<1>(N) * std::get<2>(N), t), M_(std::get<0>(N)), N_(std::get<1>(N)),
        P_(std::get<2>(N)) {}

  // forward construction
  Array(int M, int N, int P, T *d) : Impl::Data_<T>(N * M * P, d), M_(M), N_(N), P_(P) {}
  Array(const std::string& name, int M, int N, int P, T *d) : name_(name), Impl::Data_<T>(N * M * P, d), M_(M), N_(N), P_(P) {}
  Array(std::array<int, 3> N, T *d)
      : Impl::Data_<T>(std::get<0>(N) * std::get<1>(N) * std::get<2>(N), d), M_(std::get<0>(N)), N_(std::get<1>(N)),
        P_(std::get<2>(N)) {}
  Array(const std::string& name, std::array<int, 3> N, T *d)
      : name_(name), Impl::Data_<T>(std::get<0>(N) * std::get<1>(N) * std::get<2>(N), d), M_(std::get<0>(N)), N_(std::get<1>(N)),
        P_(std::get<2>(N)) {}


  // forward construction
  Array(const Array<T, 3> &other) = default;

  T &operator()(int i, int j, int k) {
    assert(0 <= i && i < M_ && 0 <= j && j < N_ && 0 <= k && k < P_);
    return d_[k + P_ * (j + i * N_)];
  }
  const T &operator()(int i, int j, int k) const {
    assert(0 <= i && i < M_ && 0 <= j && j < N_ && 0 <= k && k < P_);
    return d_[k + P_ * (j + i * N_)];
  }

  Array<T, 2> operator[](int i) {
    assert(0 <= i && i < M_);
    return Array<T, 2>(name_, N_, P_, &d_[i * N_ * P_]);
  }
  const Array<T, 2> operator[](int i) const {
    assert(0 <= i && i < M_);
    return Array<T, 2>(name_, N_, P_, &d_[i * N_ * P_]);
  }

  // resize accessor - destructive
  Array<T, 3>& resize(int M, int N, int P) {
    Impl::Data_<T>::Data_resize(M*N*P);
    M_ = M;
    N_ = N;
    P_ = P;
    return *this;
  }

  // shape
  int extent(int d) const {
    assert(d < 3 && "Array::extent requested for dimension greater than is stored.");
    return shape()[d];
  }

  std::array<int, 3> shape() const { return {M_, N_, P_}; }

  // get variable name
  std::string label() const { return {name_}; }

  typedef T value_type;

protected:
  int M_, N_, P_;
  using Impl::Data_<T>::d_;
  const std::string name_ = ("");
};

// 4D specialization
template <typename T> class Array<T, 4> : public Impl::Data_<T> {
  static const size_t dim = 4;

public:
  // forward construction
  Array(int M, int N, int P, int Q) : Impl::Data_<T>(N * M * P * Q), M_(M), N_(N), P_(P), Q_(Q) {}
  Array(std::array<int, 4> N)
      : Impl::Data_<T>(std::get<0>(N) * std::get<1>(N) * std::get<2>(N) * std::get<3>(N)), M_(std::get<0>(N)),
        N_(std::get<1>(N)), P_(std::get<2>(N)), Q_(std::get<3>(N)) {}

  // forward construction
  Array(int M, int N, int P, int Q, T t) : Impl::Data_<T>(N * M * P * Q, t), M_(M), N_(N), P_(P), Q_(Q) {}
  Array(std::array<int, 4> N, T t)
      : Impl::Data_<T>(std::get<0>(N) * std::get<1>(N) * std::get<2>(N) * std::get<3>(N), t), M_(std::get<0>(N)),
        N_(std::get<1>(N)), P_(std::get<2>(N)), Q_(std::get<3>(N)) {}

  // forward construction
  Array(int M, int N, int P, int Q, T *d) : Impl::Data_<T>(N * M * P * Q, d), M_(M), N_(N), P_(P), Q_(Q) {}
  Array(std::array<int, 4> N, T *d)
      : Impl::Data_<T>(std::get<0>(N) * std::get<1>(N) * std::get<2>(N) * std::get<3>(N), d), M_(std::get<0>(N)),
        N_(std::get<1>(N)), P_(std::get<2>(N)), Q_(std::get<3>(N)) {}

  // forward construction
  Array(const Array<T, 4> &other) = default;

  T &operator()(int i, int j, int k, int l) {
    assert(0 <= i && i < M_ && 0 <= j && j < N_ && 0 <= k && k < P_ && 0 <= l && l < Q_);
    return d_[l + Q_ * (k + P_ * (j + i * N_))];
  }
  const T &operator()(int i, int j, int k, int l) const {
    assert(0 <= i && i < M_ && 0 <= j && j < N_ && 0 <= k && k < P_ && 0 <= l && l < Q_);
    return d_[l + Q_ * (k + P_ * (j + i * N_))];
  }

  Array<T, 3> operator[](int i) {
    assert(0 <= i && i < M_);
    return Array<T, 3>(N_, P_, Q_, &d_[i * N_ * P_ * Q_]);
  }
  const Array<T, 3> operator[](int i) const {
    assert(0 <= i && i < M_);
    return Array<T, 3>(N_, P_, Q_, &d_[i * N_ * P_ * Q_]);
  }

  // resize accessor - destructive
  Array<T, 4>& resize(int M, int N, int P, int Q) {
    Impl::Data_<T>::Data_resize(M*N*P*Q);
    M_ = M;
    N_ = N;
    P_ = P;
    Q_ = Q;
    return *this;
  }

  // shape
  int extent(int d) const {
    assert(d < 4 && "Array::extent requested for dimension greater than is stored.");
    return shape()[d];
  }

  std::array<int, 4> shape() const { return {M_, N_, P_, Q_}; }

protected:
  int M_, N_, P_, Q_;
  using Impl::Data_<T>::d_;
};

//
// Construct a non-owning Array of a different shape but with the same data as
// an old shape.
//
// NOTE, this does no transposing!
//
template <size_t D1, size_t D2, typename T>
Array<T, D2> reshape(Array<T, D1> &arr_in, const std::array<int, D2> &new_shape) {
  int new_length = std::accumulate(new_shape.begin(), new_shape.end(), 1, std::multiplies<int>());
  assert(new_length == arr_in.size() && "Invalid Array reshape");

  // construct and return
  return Array<T, D2>(new_shape, (T *)arr_in.begin());
}

//
// Copies from an Array into an Array-like object.
//
template <typename T, typename Array_type> void deep_copy(Array_type &arr, const Array<T, 1> &arr_in) {
  assert(arr.extent(0) == arr_in.extent(0));
  for (int i = 0; i != arr_in.extent(0); ++i) {
    arr(i) = arr_in(i);
  }
}

template <typename T, typename Array_type> void deep_copy(Array_type &arr, const Array<T, 2> &arr_in) {
  assert(arr.extent(0) == arr_in.extent(0));
  assert(arr.extent(1) == arr_in.extent(1));
  for (int i = 0; i != arr_in.extent(0); ++i) {
    for (int j = 0; j != arr_in.extent(1); ++j) {
      arr(i, j) = arr_in(i, j);
    }
  }
}

template <typename T, typename Array_type> void deep_copy(Array_type &arr, const Array<T, 3> &arr_in) {
  assert(arr.extent(0) == arr_in.extent(0));
  assert(arr.extent(1) == arr_in.extent(1));
  assert(arr.extent(2) == arr_in.extent(2));
  for (int i = 0; i != arr_in.extent(0); ++i) {
    for (int j = 0; j != arr_in.extent(1); ++j) {
      for (int k = 0; k != arr_in.extent(2); ++k) {
        arr(i, j, k) = arr_in(i, j, k);
      }
    }
  }
}

template <typename T, typename Array_type> void deep_copy(Array_type &arr, const Array<T, 4> &arr_in) {
  assert(arr.extent(0) == arr_in.extent(0));
  assert(arr.extent(1) == arr_in.extent(1));
  assert(arr.extent(2) == arr_in.extent(2));
  assert(arr.extent(3) == arr_in.extent(3));
  for (int i = 0; i != arr_in.extent(0); ++i) {
    for (int j = 0; j != arr_in.extent(1); ++j) {
      for (int k = 0; k != arr_in.extent(2); ++k) {
        for (int l = 0; l != arr_in.extent(3); ++l) {
          arr(i, j, k, l) = arr_in(i, j, k, l);
        }
      }
    }
  }
}

template <typename T> void deep_copy(Array<T, 1> &arr, T val) {
  for (int i = 0; i != arr.extent(0); ++i)
    arr(i) = val;
}
template <typename T> void deep_copy(Array<T, 2> &arr, T val) {
  for (int i = 0; i != arr.extent(0); ++i)
    for (int j = 0; j != arr.extent(1); ++j)
      arr(i, j) = val;
}
template <typename T> void deep_copy(Array<T, 3> &arr, T val) {
  for (int i = 0; i != arr.extent(0); ++i)
    for (int j = 0; j != arr.extent(1); ++j)
      for (int k = 0; k != arr.extent(2); ++k)
        arr(i, j, k) = val;
}
template <typename T> void deep_copy(Array<T, 4> &arr, T val) {
  for (int i = 0; i != arr.extent(0); ++i)
    for (int j = 0; j != arr.extent(1); ++j)
      for (int k = 0; k != arr.extent(2); ++k)
        for (int l = 0; l != arr.extent(3); ++l)
          arr(i, j, k, l) = val;
}

} // namespace ELM

#endif
