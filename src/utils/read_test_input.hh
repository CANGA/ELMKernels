/*! \file read_test_input.h
\brief Contains utility class to read single-column input data from ELM modules
*/

#include <iostream>
#include <iterator>
#include <assert.h>

#include "array.hh"
#include <iomanip>


namespace ELM {
namespace IO {

template <typename T>
bool IsAlmostEqual(const T &a, const T &b, const double rel_tol) {
  if (a == b) return true;
  double abs_tol = 1e-20;
  auto diff = std::abs(a - b);
  auto maxreldiff = std::max(std::abs(a),std::abs(b)) * rel_tol;
  return (diff <= maxreldiff || diff <= abs_tol);
}

/*! Class to read single-column inut data from ELM modules */
class ELMtestinput {

public:
  ELMtestinput(const std::string& filename) : filename_(filename) { this->filestring_ = readInputFiletoString(); }

/*! Get state for time nstep and store in class variable state_ */
  std::string getState (const int nstep);

/*! Print state for current nstep */
  void printState ();

/*! Print entire input file */
  void printFile ();

/*! Parse state_ and assign data to arr */
  template <class Array_t>
  void parseState(Array_t arr) {
    std::string line_str, namefromfile;
    std::istringstream state_ss(this->state_);
    const std::string namefromarr = arr.label();
    while (std::getline(state_ss, line_str)) {
      std::istringstream line_ss(line_str);
      line_ss >> namefromfile;
      if (namefromfile == namefromarr) {
        const std::size_t pos = line_ss.tellg();
        const auto len = std::distance(std::istream_iterator<std::string>(line_ss), std::istream_iterator<std::string>());
        if (len != arr.size()) {
          std::string err = "INPUT ERROR: Array length (" + std::to_string(arr.size()) + ") != input data length (" + 
          std::to_string(len) + ") for variable " + namefromarr;
          throw std::runtime_error(err);
        }
        line_ss.clear();
        line_ss.seekg(pos);
        for (auto& val : arr) {
          line_ss >> val;
        }
        return;
      }
    }
    std::string err = "INPUT ERROR: Can't find variable " + namefromarr + " in NSTEP " + std::to_string(this->nstep_);
    throw std::runtime_error(err);
  }

  template <class Array_t>
  void compareOutput(Array_t arr, const double rel_tol=1e-15) {
    using ArrType = typename Array_t::value_type;
    auto filearr = Array<ArrType, 1>(arr.label(), arr.size());
    parseState(filearr);
    bool same = true;
    std::vector<std::tuple<int, ArrType, ArrType>> mismatch;
    for (std::size_t i = 0; i < arr.size(); ++i) {
      if (!IsAlmostEqual(arr.data()[i], filearr.data()[i], rel_tol)) {
        same = false;
        mismatch.push_back(std::tuple<int, ArrType, ArrType>(i, arr.data()[i], filearr.data()[i]));
      }
    }
    std::cout << std::boolalpha << arr.label() << " from NSTEP " << this->nstep_ << " passes: " << same << std::endl;
    if (!same) {
      for (const auto [i,a,f] : mismatch) {
        std::cout << std::setprecision (15) << "    [i, ELM Kernels, file]        " << i << "  " << a << "  " << f << std::endl;
      }
    }
  }


private:

/*! Read entire input file and store in class variable filestring_ */
  std::string readInputFiletoString();

  std::string filename_;
  std::string filestring_;
  std::string state_;
  int nstep_ = 0;
};

} // namespace IO
} // namespace ELM

