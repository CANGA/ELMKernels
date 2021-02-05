#pragma once

namespace ELM {

struct LandType {

  LandType() : ltype(0), ctype(0), vtype(0), urbpoi(false), lakpoi(false) {}

  // land unit, urban unit, vegetation unit
  int ltype, ctype, vtype;
  // urban point, lake point
  bool urbpoi, lakpoi;
};

} // namespace ELM
