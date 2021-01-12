#ifndef LANDTYPE_H_
#define LANDTYPE_H_

struct LandType {
  // land unit, urban unit, vegetation unit
  int ltype, ctype, vtype;
  // urban point, lake point
  bool urbpoi, lakpoi;
};

#endif
