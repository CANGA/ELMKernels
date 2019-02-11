//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_UTILS_HH_
#define ELM_KERNEL_TEST_UTILS_HH_

namespace ELM {
namespace Utils {


enum struct Ordering { C, FORTRAN };

template <typename T=double, Ordering O=Ordering::C>
class Matrix {

 public:
  Matrix(size_t nrows, size_t ncols) :
      nrows_(nrows),
      ncols_(ncols),
      d_(nrows*ncols) {}

  inline
  const T& operator()(size_t i, size_t j) const {
    return d_[ncols_*i + j];
  }

  inline
  T& operator()(size_t i, size_t j) {
    return d_[ncols_*i + j];
  }

  void operator=(T val) {
    for (size_t lcv=0; lcv!=nrows_*ncols_; ++lcv) {
      d_[lcv] = val;
    }
  }

  typename std::vector<T>::const_iterator begin() const {
    return d_.begin();
  }

  typename std::vector<T>::const_iterator end() const {
    return d_.end();
  }
  
  /*
  inline
  T* operator[](size_t i) {
    return &d_[i*ncols];
  }
  */
 private:
  size_t nrows_, ncols_;
  std::vector<T> d_;

};


template <typename T>
class Matrix<T,Ordering::FORTRAN> {

 public:
  Matrix(size_t nrows, size_t ncols) :
      nrows_(nrows),
      ncols_(ncols),
      d_(nrows*ncols) {}

  inline
  const T& operator()(size_t i, size_t j) const {
    return d_[nrows_*j + i];
  }

  inline
  T& operator()(size_t i, size_t j) {
    return d_[nrows_*j + i];
  }

  /*
  inline
  T* operator[](size_t i) {
    return &d_[i*ncols];
  }
  */
  
  void operator=(T val) {
    for (size_t lcv=0; lcv!=nrows_*ncols_; ++lcv) {
      d_[lcv] = val;
    }
  }


  typename std::vector<T>::const_iterator begin() const {
    return d_.begin();
  }

  typename std::vector<T>::const_iterator end() const {
    return d_.end();
  }
  
 private:
  size_t nrows_, ncols_;
  std::vector<T> d_;

};



template <typename T=double, Ordering O=Ordering::C>
class TensorRank3 {

 public:
  TensorRank3(size_t n_first, size_t n_second, size_t n_third) :
      n_first_(n_first),
      n_second_(n_second),
      n_third_(n_third),
      d_(n_first*n_second*n_third) {}

  inline
  const T& operator()(size_t i, size_t j, size_t k) const {
    return d_[n_third_ * (n_second_*i + j) + k];
  }

  inline
  T& operator()(size_t i, size_t j, size_t k) {
    return d_[n_third_ * (n_second_*i + j) + k];
  }

  /*
  inline
  T* operator[](size_t i) {
    return &d_[i*n_second_*n_third_];
  }
  */
  
  void operator=(T val) {
    for (size_t lcv=0; lcv!=n_first_*n_second_*n_third_; ++lcv) {
      d_[lcv] = val;
    }
  }


  typename std::vector<T>::const_iterator begin() const {
    return d_.begin();
  }

  typename std::vector<T>::const_iterator end() const {
    return d_.end();
  }
  
 private:
  size_t n_first_, n_second_, n_third_;
  std::vector<T> d_;

};


template <typename T>
class TensorRank3<T,Ordering::FORTRAN> {

 public:
  TensorRank3(size_t n_first, size_t n_second, size_t n_third) :
      n_first_(n_first),
      n_second_(n_second),
      n_third_(n_third),
      d_(n_first*n_second*n_third) {}

  inline
  const T& operator()(size_t i, size_t j, size_t k) const {
    return d_[n_first_*(n_second_*k + j) + i];
  }

  inline
  T& operator()(size_t i, size_t j, size_t k) {
    return d_[n_first_*(n_second_*k + j) + i];
  }

  /*
  inline
  T* operator[](size_t k) {
    return &d_[k*n_first_*n_second_];
  }
  */

  void operator=(T val) {
    for (size_t lcv=0; lcv!=n_first_*n_second_*n_third_; ++lcv) {
      d_[lcv] = val;
    }
  }
  
  typename std::vector<T>::const_iterator begin() const {
    return d_.begin();
  }

  typename std::vector<T>::const_iterator end() const {
    return d_.end();
  }

 private:
  size_t n_first_, n_second_, n_third_;
  std::vector<T> d_;

};

  
} // namespace Utils
} // namespace ELM


#endif
