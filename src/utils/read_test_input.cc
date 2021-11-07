#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <assert.h>

#include "array.hh"
#include "read_test_input.hh"

namespace ELM {
namespace IO {

  std::string ELMtestinput::getState (const int nstep) {
    this->nstep_ = nstep;
    const std::string startstr = "NSTEP " + std::to_string(nstep) + "\n";
    const std::string endstr = "!!! " + std::to_string(nstep) + "\n";
    const std::size_t nstepstart = filestring_.find(startstr);
    const std::size_t nstepend = filestring_.find(endstr);
    return this->state_ = this->filestring_.substr(nstepstart, nstepend - nstepstart);
  }

  void ELMtestinput::printState () {
     std::cout << "State for NSTEP " << this->nstep_ << ": " << this->state_ << std::endl;
  }

  void ELMtestinput::printFile () {
     std::cout << this->filename_ << " Input file contains: " << this->filestring_ << std::endl;
  }

  std::string ELMtestinput::readInputFiletoString() {
    std::string in_string;
    std::ifstream in_file(this->filename_);
    if (!in_file) {std::string err = "INPUT ERROR: Can't open input file " + this->filename_; throw std::runtime_error(err);}
    in_file.seekg(0, std::ios::end);
    in_string.resize(in_file.tellg());
    in_file.seekg(std::ios::beg);
    in_file.read(&in_string[0], in_string.size());
    in_file.close();
    return this->filestring_ = in_string;
  }

} // namespace IO
} // namespace ELM

