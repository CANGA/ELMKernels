#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

class ELMtestinput {

public:
  ELMtestinput(const std::string filename) : filename_(filename) { filestring_ = readInputFiletoString(); }


  std::string getState (const int nstep) {
    this->nstep_ = nstep;
    std::string startstr = "NSTEP " + std::to_string(nstep) + "\n";
    std::string endstr = "!!! " + std::to_string(nstep) + "\n";
    int nstepstart = filestring_.find(startstr);
    int nstepend = filestring_.find(endstr);
    return state_ = filestring_.substr(nstepstart, nstepend - nstepstart);
  }

  template <class Array_t>
  int parseState(Array_t& arr) {
    std::string text_read, varname, data;
    std::istringstream state(state_);
    while (std::getline(state, text_read))
    {
      std::istringstream ss(text_read);
      ss >> varname;
      if (varname == arr.getname()) {
        for (auto& d : arr) {
          ss >> d;
        }
        return 0;
      }
    }
    return 1;
  }

  void printState () {
     std::cout << "State for NSTEP " << nstep_ << ": " << state_ << std::endl;
  }

  void printFile () {
     std::cout << filename_ << " Contains: " << filestring_ << std::endl;
  }

  bool checkStateNstep () {
    std::string readnstep;
    std::istringstream state(state_);
    std::getline(state, readnstep);
    std::string compstr = "NSTEP " + std::to_string(nstep_);
    return (readnstep == compstr);
  }

private:

  std::string readInputFiletoString() {
    std::ifstream in_file(filename_);
    std::string in_string;
    if (!in_file) {std::string err = "INPUT ERROR: Can't open input file " + filename_; throw std::runtime_error(err);}
    in_file.seekg(0, std::ios::end);
    in_string.resize(in_file.tellg());
    in_file.seekg(std::ios::beg);
    in_file.read(&in_string[0], in_string.size());
    in_file.close();
    return filestring_ = in_string;
  }

  std::string filename_;
  std::string filestring_;
  std::string state_;
  int nstep_ = 0;
};
