//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_UTILS_HH_
#define ELM_KERNEL_TEST_UTILS_HH_

#include <iostream>
#include <memory>
#include <type_traits>

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
  
 protected:
  //
  // Constructors are all protected -- do not use this class directly!
  //
  // construct from a total length
  Data_(std::size_t len)
      : len_(len),
        do_(std::shared_ptr<T>(new T[len], std::default_delete<T[]>()))
  {
    d_ = do_.get();
  }

  // construct and initialize
  Data_(std::size_t len, T d)
      : Data_(len)
  {
    *this = d;
  }

  // construct non-owning view
  Data_(std::size_t len, T* d)
      : len_(len),
        d_(d) {}

  // copy construct  
  Data_(const Data_& other) = default;

 protected:
  // global length
  const std::size_t len_;

  // owning data
  std::shared_ptr<T> do_;

  // non-owning data
  T* d_;
};
  

// templated class for multi-dimensional array
template<int D, typename T>
class Array : public Data_<T> {};

// 1D specialization
template<typename T>
class Array<1,T> : public Data_<T> {
 public:
  // forward construction
  Array(std::size_t N)
      : Data_<T>(N) {}

  // forward construction
  Array(std::size_t N, T t)
      : Data_<T>(N, t) {}
  
  // forward construction
  Array(std::size_t N, T* d) :
      Data_<T>(N, d)
  {}

  // forward construction
  Array(const Array<1,T>& other) = default;

  // accessors
  T& operator()(const std::size_t i) { assert(0 <= i && i < len_); return d_[i]; }
  const T& operator()(const std::size_t i) const { assert(0 <= i && i < len_); return d_[i]; }

  T& operator[](const std::size_t i) { assert(0 <= i && i < len_); return d_[i]; }
  const T& operator[](std::size_t i) const { assert(0 <= i && i < len_); return d_[i]; }

 protected:
  using Data_<T>::len_;
  using Data_<T>::d_;
};


// 2D specialization
template<typename T>
class Array<2,T> : public Data_<T> {
 public:
  // forward construction
  Array(std::size_t M, std::size_t N)
      : Data_<T>(N*M),
        M_(M),
        N_(N)
  {}

  // forward construction
  Array(std::size_t M, std::size_t N, T t)
      : Data_<T>(N*M, t),
        M_(M),
        N_(N)
  {}
  
  // forward construction
  Array(std::size_t M, std::size_t N, T* d)
      : Data_<T>(N*M, d),
        M_(M),
        N_(N)
  {}

  // forward construction
  Array(const Array<2,T>& other) = default;

  T& operator()(std::size_t i, std::size_t j) { assert(0 <= i && i < M_ && 0 <= j && j < N_); return d_[j+i*N_]; }
  const T& operator()(std::size_t i, std::size_t j) const { assert(0 <= i && i < M_ && 0 <= j && j < N_); return d_[j+i*N_]; }

  Array<1,T> operator[](std::size_t i) { assert(0 <= i && i < M_); return Array<1,T>(N_, &d_[i*N_]); }
  const Array<1,T> operator[](std::size_t i) const { assert(0 <= i && i < M_); return Array<1,T>(N_, &d_[i*N_]); }
  
 protected:
  std::size_t M_, N_;
  using Data_<T>::d_;

};


// 3D specialization
template<typename T>
class Array<3,T> : public Data_<T> {
 public:
  // forward construction
  Array(std::size_t M, std::size_t N, std::size_t P)
      : Data_<T>(N*M*P),
        M_(M),
        N_(N),
        P_(P)
  {}

  // forward construction
  Array(std::size_t M, std::size_t N, std::size_t P, T t)
      : Data_<T>(N*M*P, t),
        M_(M),
        N_(N),
        P_(P)
  {}
  
  // forward construction
  Array(std::size_t M, std::size_t N, std::size_t P, T* d)
      : Data_<T>(N*M*P, d),
        M_(M),
        N_(N),
        P_(P)
  {}

  // forward construction
  Array(const Array<3,T>& other) = default;

  T& operator()(std::size_t i, std::size_t j, std::size_t k) {
    assert(0 <= i && i < M_ &&
           0 <= j && j < N_ &&
           0 <= k && k < P_);
    return d_[k + P_*(j+i*N_)];
  }
  const T& operator()(std::size_t i, std::size_t j, std::size_t k) const {
    assert(0 <= i && i < M_ &&
           0 <= j && j < N_ &&
           0 <= k && k < P_);
    return d_[k + P_*(j+i*N_)];
  }

  Array<2,T> operator[](std::size_t i) {
    assert(0 <= i && i < M_);
    return Array<2,T>(N_, P_, &d_[i*N_*P_]);
  }
  const Array<2,T> operator[](std::size_t i) const {
    assert(0 <= i && i < M_);
    return Array<2,T>(N_, P_, &d_[i*N_*P_]);
  }
  
 protected:
  std::size_t M_, N_, P_;
  using Data_<T>::d_;

};



// 4D specialization
template<typename T>
class Array<4,T> : public Data_<T> {
 public:
  // forward construction
  Array(std::size_t M, std::size_t N, std::size_t P, std::size_t Q)
      : Data_<T>(N*M*P*Q),
        M_(M),
        N_(N),
        P_(P),
        Q_(Q)
  {}

  // forward construction
  Array(std::size_t M, std::size_t N, std::size_t P, std::size_t Q, T t)
      : Data_<T>(N*M*P*Q, t),
        M_(M),
        N_(N),
        P_(P),
        Q_(Q)
  {}
  
  // forward construction
  Array(std::size_t M, std::size_t N, std::size_t P, std::size_t Q, T* d)
      : Data_<T>(N*M*P*Q, d),
        M_(M),
        N_(N),
        P_(P),
        Q_(Q)
  {}

  // forward construction
  Array(const Array<4,T>& other) = default;

  T& operator()(std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    assert(0 <= i && i < M_ &&
           0 <= j && j < N_ &&
           0 <= k && k < P_ &&
           0 <= l && l < Q_);
    return d_[l + Q_*(k + P_*(j+i*N_))];
  }
  const T& operator()(std::size_t i, std::size_t j, std::size_t k) const {
    assert(0 <= i && i < M_ &&
           0 <= j && j < N_ &&
           0 <= k && k < P_ &&
           0 <= l && l < Q_);
    return d_[l + Q_*(k + P_*(j+i*N_))];
  }

  Array<3,T> operator[](std::size_t i) {
    assert(0 <= i && i < M_);
    return Array<3,T>(N_, P_, Q_, &d_[i*N_*P_*Q_]);
  }
  const Array<3,T> operator[](std::size_t i) const {
    assert(0 <= i && i < M_);
    return Array<3,T>(N_, P_, Q_, &d_[i*N_*P_*Q_]);
  }

  
 protected:
  std::size_t M_, N_, P_, Q_;
  using Data_<T>::d_;

};


  
} // namespace Utils
} // namespace ELM


#endif
