//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_ARRAY_HH_
#define ELM_KERNEL_TEST_ARRAY_HH_

#include <tuple>
#include <algorithm>
#include <array>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <assert.h>
#include <memory>

namespace ELM {
namespace Utils {


template<typename T>
class Data_ {
 public:
  // assigment operator
  void operator=(T t) {
    for (auto& d : *this) {
      d = t;
    }
  }

  // iterators
  T const * begin() const { return &(d_[0]); }
  T const * end() const { return &(d_[len_]); }
  T * begin() { return &(d_[0]); }
  T * end() { return &(d_[len_]); }
  size_t size() const { return len_; }
  
 protected:
  //
  // Constructors are all protected -- do not use this class directly!
  //
  // construct from a total length
  Data_(size_t len)
      : len_(len),
        do_(std::shared_ptr<T>(new T[len], std::default_delete<T[]>()))
  {
    d_ = do_.get();
  }

  // construct and initialize
  Data_(size_t len, T d)
      : Data_(len)
  {
    *this = d;
  }

  // construct non-owning view
  Data_(size_t len, T* d)
      : len_(len),
        d_(d) {}

  // copy construct  
  Data_(const Data_& other) = default;

 protected:
  // global length
  const size_t len_;

  // owning data
  std::shared_ptr<T> do_;

  // non-owning data
  T* d_;
};
  

// templated class for multi-dimensional array
template<typename T, size_t D>
class Array : public Data_<T> {};

// 1D specialization
template<typename T>
class Array<T,1> : public Data_<T> {
  static const size_t dim = 1;

 public:
  // forward construction
  Array(size_t N)
      : Data_<T>(N) {}

  // forward construction
  Array(size_t N, T t)
      : Data_<T>(N, t) {}
  
  // forward construction
  Array(size_t N, T* d) :
      Data_<T>(N, d)
  {}

  // forward construction
  Array(const Array<T,1>& other) = default;

  // accessors
  T& operator()(const size_t i) { assert(0 <= i && i < len_); return d_[i]; }
  const T& operator()(const size_t i) const { assert(0 <= i && i < len_); return d_[i]; }

  T& operator[](const size_t i) { assert(0 <= i && i < len_); return d_[i]; }
  const T& operator[](size_t i) const { assert(0 <= i && i < len_); return d_[i]; }

  // shape
  size_t extent(size_t d) const {
    assert(d < 1 && "Array::extent requested for dimension greater than is stored.");
    return shape()[d];
  }
  
  std::array<size_t,1> shape() const {
    return { len_ };
  }
  
 protected:
  using Data_<T>::len_;
  using Data_<T>::d_;
};


// 2D specialization
template<typename T>
class Array<T,2> : public Data_<T> {
  static const size_t dim = 2;

 public:
  // forward construction
  Array(size_t M, size_t N)
      : Data_<T>(N*M),
        M_(M),
        N_(N)
  {}

  // forward construction
  Array(size_t M, size_t N, T t)
      : Data_<T>(N*M, t),
        M_(M),
        N_(N)
  {}
  
  // forward construction
  Array(size_t M, size_t N, T* d)
      : Data_<T>(N*M, d),
        M_(M),
        N_(N)
  {}

  // forward construction
  Array(const Array<T,2>& other) = default;

  T& operator()(size_t i, size_t j) { assert(0 <= i && i < M_ && 0 <= j && j < N_); return d_[j+i*N_]; }
  const T& operator()(size_t i, size_t j) const { assert(0 <= i && i < M_ && 0 <= j && j < N_); return d_[j+i*N_]; }

  Array<T,1> operator[](size_t i) { assert(0 <= i && i < M_); return Array<T,1>(N_, &d_[i*N_]); }
  const Array<T,1> operator[](size_t i) const { assert(0 <= i && i < M_); return Array<T,1>(N_, &d_[i*N_]); }

  // shape
  size_t extent(size_t d) const {
    assert(d < 2 && "Array::extent requested for dimension greater than is stored.");
    return shape()[d];
  }
  
  std::array<size_t,2> shape() const {
    return { M_, N_ };
  }

  
 protected:
  size_t M_, N_;
  using Data_<T>::d_;

};


// 3D specialization
template<typename T>
class Array<T,3> : public Data_<T> {
  static const size_t dim = 3;
  
 public:
  // forward construction
  Array(size_t M, size_t N, size_t P)
      : Data_<T>(N*M*P),
        M_(M),
        N_(N),
        P_(P)
  {}

  // forward construction
  Array(size_t M, size_t N, size_t P, T t)
      : Data_<T>(N*M*P, t),
        M_(M),
        N_(N),
        P_(P)
  {}
  
  // forward construction
  Array(size_t M, size_t N, size_t P, T* d)
      : Data_<T>(N*M*P, d),
        M_(M),
        N_(N),
        P_(P)
  {}

  // forward construction
  Array(const Array<T,3>& other) = default;

  T& operator()(size_t i, size_t j, size_t k) {
    assert(0 <= i && i < M_ &&
           0 <= j && j < N_ &&
           0 <= k && k < P_);
    return d_[k + P_*(j+i*N_)];
  }
  const T& operator()(size_t i, size_t j, size_t k) const {
    assert(0 <= i && i < M_ &&
           0 <= j && j < N_ &&
           0 <= k && k < P_);
    return d_[k + P_*(j+i*N_)];
  }

  Array<T,2> operator[](size_t i) {
    assert(0 <= i && i < M_);
    return Array<T,2>(N_, P_, &d_[i*N_*P_]);
  }
  const Array<T,2> operator[](size_t i) const {
    assert(0 <= i && i < M_);
    return Array<T,2>(N_, P_, &d_[i*N_*P_]);
  }

  // shape
  size_t extent(size_t d) const {
    assert(d < 3 && "Array::extent requested for dimension greater than is stored.");
    return shape()[d];
  }
  
  std::array<size_t,3> shape() const {
    return { M_, N_, P_ };
  }
  
 protected:
  size_t M_, N_, P_;
  using Data_<T>::d_;

};



// 4D specialization
template<typename T>
class Array<T,4> : public Data_<T> {
  static const size_t dim = 4;
 public:
  // forward construction
  Array(size_t M, size_t N, size_t P, size_t Q)
      : Data_<T>(N*M*P*Q),
        M_(M),
        N_(N),
        P_(P),
        Q_(Q)
  {}

  // forward construction
  Array(size_t M, size_t N, size_t P, size_t Q, T t)
      : Data_<T>(N*M*P*Q, t),
        M_(M),
        N_(N),
        P_(P),
        Q_(Q)
  {}
  
  // forward construction
  Array(size_t M, size_t N, size_t P, size_t Q, T* d)
      : Data_<T>(N*M*P*Q, d),
        M_(M),
        N_(N),
        P_(P),
        Q_(Q)
  {}

  // forward construction
  Array(const Array<T,4>& other) = default;

  T& operator()(size_t i, size_t j, size_t k, size_t l) {
    assert(0 <= i && i < M_ &&
           0 <= j && j < N_ &&
           0 <= k && k < P_ &&
           0 <= l && l < Q_);
    return d_[l + Q_*(k + P_*(j+i*N_))];
  }
  const T& operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(0 <= i && i < M_ &&
           0 <= j && j < N_ &&
           0 <= k && k < P_ &&
           0 <= l && l < Q_);
    return d_[l + Q_*(k + P_*(j+i*N_))];
  }

  Array<T,3> operator[](size_t i) {
    assert(0 <= i && i < M_);
    return Array<T,3>(N_, P_, Q_, &d_[i*N_*P_*Q_]);
  }
  const Array<T,3> operator[](size_t i) const {
    assert(0 <= i && i < M_);
    return Array<T,3>(N_, P_, Q_, &d_[i*N_*P_*Q_]);
  }

  // shape
  size_t extent(size_t d) const {
    assert(d < 4 && "Array::extent requested for dimension greater than is stored.");
    return shape()[d];
  }
  
  std::array<size_t,4> shape() const {
    return { M_, N_, P_, Q_ };
  }
  
 protected:
  size_t M_, N_, P_, Q_;
  using Data_<T>::d_;

};


//
// Construct a non-owning Array of a different shape but with the same data as
// an old shape.
//
// NOTE, this does no transposing!
//
template<size_t D1, size_t D2, typename T>
Array<T,D2> reshape(Array<T,D1>& arr_in, const std::array<size_t, D2>& new_shape) {
  size_t new_length = std::accumulate(new_shape.begin(), new_shape.end(), 1, std::multiplies<size_t>());
  if (new_length != arr_in.size()) {
    //std::stringstream err;
    std::cout << "Invalid Array reshape, cannot reshape object of size: " << arr_in.size() << " into array of size: " << new_length;
  }

  // make a tuple of the new shape plus the pointer to data.
  auto data_tuple = std::make_tuple( (T*) arr_in.begin() );
  auto constructor_args = std::tuple_cat(new_shape, data_tuple);

  // construct and return
  return std::make_from_tuple<Array<T,D2>>(std::move(constructor_args));
}


//
// Copies from an Array into an Array-like object.
//
template<typename T, typename Array_type>
void
deep_copy(Array_type& arr, const Array<T,1>& arr_in)
{
  assert(arr.extent(0) == arr_in.extent(0));
  for (size_t i=0; i!=arr_in.extent(0); ++i) {
    arr(i) = arr_in(i);
  }
}

template<typename T, typename Array_type>
void
deep_copy(Array_type& arr, const Array<T,2>& arr_in)
{
  assert(arr.extent(0) == arr_in.extent(0));
  assert(arr.extent(1) == arr_in.extent(1));
  for (size_t i=0; i!=arr_in.extent(0); ++i) {
    for (size_t j=0; j!=arr_in.extent(1); ++j) {
      arr(i,j) = arr_in(i,j);
    }
  }
}

template<typename T, typename Array_type>
void
deep_copy(Array_type& arr, const Array<T,3>& arr_in)
{
  assert(arr.extent(0) == arr_in.extent(0));
  assert(arr.extent(1) == arr_in.extent(1));
  assert(arr.extent(2) == arr_in.extent(2));
  for (size_t i=0; i!=arr_in.extent(0); ++i) {
    for (size_t j=0; j!=arr_in.extent(1); ++j) {
      for (size_t k=0; k!=arr_in.extent(2); ++k) {
        arr(i,j,k) = arr_in(i,j,k);
      }
    }
  }
}

template<typename T, typename Array_type>
void
deep_copy(Array_type& arr, const Array<T,4>& arr_in)
{
  assert(arr.extent(0) == arr_in.extent(0));
  assert(arr.extent(1) == arr_in.extent(1));
  assert(arr.extent(2) == arr_in.extent(2));
  assert(arr.extent(3) == arr_in.extent(3));
  for (size_t i=0; i!=arr_in.extent(0); ++i) {
    for (size_t j=0; j!=arr_in.extent(1); ++j) {
      for (size_t k=0; k!=arr_in.extent(2); ++k) {
        for (size_t l=0; l!=arr_in.extent(3); ++l) {
          arr(i,j,k,l) = arr_in(i,j,k,l);
        }
      }
    }
  }
}


} // namespace Utils
} // namespace ELM


#endif
