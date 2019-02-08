module column_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing landunit indices and associated variables and routines.
  !
  ! Drastically simplified from the original version --etc
  !-----------------------------------------------------------------------
  use landunit_varcon, only : isturb_MIN

  implicit none
  save
  private

  !------------------------------------------------------------------
  ! Initialize column type constants
  !------------------------------------------------------------------

  ! urban column types

  integer, parameter, public :: icol_roof        = isturb_MIN*10 + 1
  integer, parameter, public :: icol_sunwall     = isturb_MIN*10 + 2
  integer, parameter, public :: icol_shadewall   = isturb_MIN*10 + 3
  integer, parameter, public :: icol_road_imperv = isturb_MIN*10 + 4
  integer, parameter, public :: icol_road_perv   = isturb_MIN*10 + 5

end module column_varcon
