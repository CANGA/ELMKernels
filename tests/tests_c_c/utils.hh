//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_UTILS_HH_
#define ELM_KERNEL_TEST_UTILS_HH_

#include <iostream>
#include <memory>
#include <type_traits>

namespace ELM {
namespace Utils {


template<typename T=double>
class Vector {
 public:
  Vector(std::size_t N) :
      N_(N),
      do_(std::shared_ptr<T>(new T[N], std::default_delete<T[]>()))
  {
    d_ = do_.get();
  }
  
  Vector(std::size_t N, T t) :
      Vector(N)
  { *this = t; }

  Vector(std::size_t N, T* d) :
      N_(N),
      do_(nullptr),
      d_(d)
  {}

  Vector(const Vector& other) = default;

  T& operator()(size_t i) { assert(0 <= i && i < N_); return d_[i]; }
  const T& operator()(size_t i) const { assert(0 <= i && i < N_); return d_[i]; }

  T& operator[](size_t i) { assert(0 <= i && i < N_); return d_[i]; }
  const T& operator[](size_t i) const { assert(0 <= i && i < N_); return d_[i]; }
  
  void operator=(T t) {
    for (size_t i=0; i!=N_; ++i) {
      d_[i] = t;
    }
  }

  double const * begin() const { return &(d_[0]); }
  double const * end() const { return &(d_[N_]); }
  std::size_t size() const { return N_; }
  
 private:
  const std::size_t N_;
  std::shared_ptr<T> do_;
  T* d_;
};



template<typename T=double>
class Matrix {
 public:
  Matrix(std::size_t M, std::size_t N) :
      M_(M),
      N_(N),
      NM_(N*M),
      d_(std::shared_ptr<T>(new T[N*M], std::default_delete<T[]>()))
  {}

  Matrix(std::size_t M, std::size_t N, T t) :
      Matrix(M,N) 
  { *this = t; }
  
  Matrix(const Matrix& other) = default;
  
  T& operator()(size_t i, size_t j) { assert(0 <= i && i < M_ && 0 <= j && j < N_); return d_.get()[j+i*N_]; }
  const T& operator()(size_t i, size_t j) const { assert(0 <= i && i < M_ && 0 <= j && j < N_); return d_.get()[j+i*N_]; }

  Vector<T> operator[](size_t i) { assert(0 <= i && i < M_); return Vector<T>(N_, &d_.get()[i*N_]); }
  const Vector<T> operator[](size_t i) const { assert(0 <= i && i < M_); return Vector<T>(N_, &d_.get()[i*N_]); }

  void operator=(T t) {
    for (std::size_t i=0; i!=NM_; ++i) {
      d_.get()[i] = t;
    }
  }

  double const * begin() const { return &d_.get()[0]; }
  double const * end() const { return &d_.get()[M_-1][N_-1] +1; }
  
 private:
  std::size_t M_, N_, NM_;
  
  std::shared_ptr<T> d_;
};
  
  
} // namespace Utils
} // namespace ELM


#endif
