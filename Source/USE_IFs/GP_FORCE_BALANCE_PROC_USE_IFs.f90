! Begin MIT license text.
! _______________________________________________________________________________________________________

! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.

      MODULE GP_FORCE_BALANCE_PROC_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE GP_FORCE_BALANCE_PROC

      USE OURTIM_Interface
      USE GET_GRID_AND_COMP_Interface
      USE MODAL_OUTPUT_WRITERS, ONLY  :  WRITE_SUBCASE_EIGENVEC_HEADER
      USE GET_GRID_NUM_COMPS_Interface
      USE GET_ARRAY_ROW_NUM_Interface
      USE GET_ELGP_Interface
      USE GET_ELEM_AGRID_BGRID_Interface
      USE EMG_Interface
      USE ELEMENT_TRANSFORMATIONS, ONLY:  ELEM_TRANSFORM_LBG
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE RESULT_COORDINATES, ONLY    :  ELMDIS, TRANSFORM_NODE_FORCES
      USE ELEMENT_RECOVERY_SUPPORT, ONLY:  CALC_ELEM_NODE_FORCES
      USE COUNTER_INIT_Interface
      USE COUNTER_PROGRESS_Interface

      END MODULE GP_FORCE_BALANCE_PROC_USE_IFs
