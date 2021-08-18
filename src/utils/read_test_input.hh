/*! \file read_test_input.h
\brief Contains utility class to read single-column inut data from ELM modules
*/

#include <iterator>
#include <assert.h>

#include "array.hh"

namespace ELM {
namespace IO {

/*! Class to read single-column inut data from ELM modules */
class ELMtestinput {

public:
  ELMtestinput(const std::string filename) : filename_(filename) { filestring_ = readInputFiletoString(); }

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
    std::istringstream state_ss(state_);
    const std::string namefromarr = arr.getname();
    while (std::getline(state_ss, line_str)) {
      std::istringstream line_ss(line_str);
      line_ss >> namefromfile;
      if (namefromfile == namefromarr) {
        const std::size_t pos = line_ss.tellg();
        const auto len = std::distance(std::istream_iterator<std::string>(line_ss), std::istream_iterator<std::string>());
        assert(len == arr.extent(0) && "INPUT ERROR: Array length != input data length");
        line_ss.clear();
        line_ss.seekg(pos);
        for (auto& val : arr) {
          line_ss >> val;
        }
        return;
      }
    }
    std::string err = "INPUT ERROR: Can't find variable " + namefromarr + " in NSTEP " + std::to_string(nstep_);
    throw std::runtime_error(err);
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

