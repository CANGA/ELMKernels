//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_UTILS_HH_
#define ELM_KERNEL_TEST_UTILS_HH_

#include <iostream>
#include <memory>
#include <type_traits>

namespace ELM {
namespace Utils {


template<typename T=double>
class VectorStatic {
 public:
  VectorStatic(std::size_t N) :
      do_(std::shared_ptr<T>(new T[N])),
      N_(N) {
    d_ = do_.get();
  }
  
  VectorStatic(std::size_t N, T t) :
      VectorStatic(N)
  { *this = t; }

  VectorStatic(std::size_t N, T* d) :
      do_(nullptr),
      d_(d),
      N_(N)
  {}

  VectorStatic(const VectorStatic& other) = default;

  T& operator()(size_t i) { assert(0 <= i && i < N); return (*d_)[i]; }
  const T& operator()(size_t i) const { assert(0 <= i && i < N); return (*d_)[i]; }

  T& operator[](size_t i) { assert(0 <= i && i < N); return (*d_)[i]; }
  const T& operator[](size_t i) const { assert(0 <= i && i < N); return (*d_)[i]; }
  
  void operator=(T t) {
    for (size_t i=0; i!=N_; ++i) {
      (*d_)[i] = t;
    }
  }

  double const * begin() const { return &((*d_)[0]); }
  double const * end() const { return &((*d_)[N]); }
  std::size_t size() const { return N_; }
  
 private:
  const std::size_t N_;
  T[] d_;
  std::shared_ptr<T[]> do_;
};



template<typename T=double>
class MatrixStatic {
 public:
  MatrixStatic(std::size_t M, std::size_t N) :
      d_(std::shared_ptr<T>(new T[N*M])),
      M_(M),
      N_(N),
      NM_(N*M)
  {}

  MatrixStatic(std::size_t M, std::size_t N, T t) :
      MatrixStatic(M,N) 
  { *this = t; }
  
  MatrixStatic(const MatrixStatic& other) = default;
  
  T& operator()(size_t i, size_t j) { assert(0 <= i && i < M_ && 0 <= j && j < N_); return (*d_)[j+i*N_]; }
  const T& operator()(size_t i, size_t j) const { assert(0 <= i && i < M_ && 0 <= j && j < N_); return (*d_)[j+i*N_]; }

  VectorStatic<T> operator[](size_t i) { assert(0 <= i && i < M_); return VectorStatic<T>(N_, &(*d_)[i*N_]); }
  const VectorStatic<COL,T> operator[](size_t i) const { assert(0 <= i && i < M_); return VectorStatic<T>(N_, &(*d_)[i*N_]); }

  void operator=(T t) {
    for (std::size_t i=0; i!=NM_; ++i) {
      (*d_)[i] = t;
    }
  }

  double const * begin() const { return &(*d_)[0]; }
  double const * end() const { return &(*d_)[M_-1][N_-1] +1; }
  
 private:
  std::shared_ptr<T[]> d_;
};
  
  
} // namespace Utils
} // namespace ELM


#endif
