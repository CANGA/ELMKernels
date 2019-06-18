#ifndef LANDUNIT_VARCON_HH_
#define LANDUNIT_VARCON_HH_

  static const int istsoil    = 1  ;
  static const int istcrop    = 2  ;
  static const int istice     = 3  ;
  static const int istice_mec = 4  ;
  static const int istdlak    = 5  ;
  static const int istwet     = 6  ;

  static const int isturb_MIN = 7  ;
  static const int isturb_tbd = 7  ;
  static const int isturb_hd  = 8  ;
  static const int isturb_md  = 9  ;
  static const int isturb_MAX = 9  ;

  static const int max_lunit  = 9  ;
                                       
  static const int landunit_name_length = 40  ;
//  character(len=landunit_name_length), public  :: landunit_names(max_lunit)  ;

  static const int numurbl = isturb_MAX - isturb_MIN + 1  ;

#endif
