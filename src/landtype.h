#pragma once

namespace ELM {

struct LandType {

  LandType() : ltype(1), ctype(0), vtype(2), urbpoi(false), lakpoi(false) {}

  // land unit, urban unit, vegetation unit
  int ltype, ctype, vtype;
  // urban point, lake point
  bool urbpoi, lakpoi;
};

} // namespace ELM
