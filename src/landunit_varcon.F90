module landunit_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing landunit indices and associated variables and routines.
  !
  ! Drastically simplifed from original --etc
  ! -----------------------------------------------------------------------
  implicit none
  save
  private
  
  !------------------------------------------------------------------
  ! Initialize landunit type constants
  !------------------------------------------------------------------

  integer, parameter, public :: istsoil    = 1  !soil         landunit type (natural vegetation)
  integer, parameter, public :: istcrop    = 2  !crop         landunit type
  integer, parameter, public :: istice     = 3  !land ice     landunit type (glacier)
  integer, parameter, public :: istice_mec = 4  !land ice (multiple elevation classes) landunit type
  integer, parameter, public :: istdlak    = 5  !deep lake    landunit type (now used for all lakes)
  integer, parameter, public :: istwet     = 6  !wetland      landunit type (swamp, marsh, etc.)

  integer, parameter, public :: isturb_MIN = 7  !minimum urban type index
  integer, parameter, public :: isturb_tbd = 7  !urban tbd    landunit type
  integer, parameter, public :: isturb_hd  = 8  !urban hd     landunit type
  integer, parameter, public :: isturb_md  = 9  !urban md     landunit type
  integer, parameter, public :: isturb_MAX = 9  !maximum urban type index

  integer, parameter, public :: max_lunit  = 9  !maximum value that lun_pp%itype can have
                                        !(i.e., largest value in the above list)

  integer, parameter, public                   :: landunit_name_length = 40  ! max length of landunit names
  character(len=landunit_name_length), public  :: landunit_names(max_lunit)  ! name of each landunit type

  ! parameters that depend on the above constants

  integer, parameter, public :: numurbl = isturb_MAX - isturb_MIN + 1   ! number of urban landunits

end module landunit_varcon
