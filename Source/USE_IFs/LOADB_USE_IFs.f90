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

      MODULE LOADB_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE LOADB

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE OUTA_HERE_Interface
      USE FFIELD_Interface
      USE FFIELD2_Interface
      USE DOF_SETS, ONLY              :  BD_ASET, BD_ASET1, BD_USET, BD_USET1, BD_SUPORT
      USE ROD_BAR_BEAM_CARDS, ONLY    :  BD_BAROR, BD_BEAMOR, BD_CBAR, BD_CROD, BD_CONROD, BD_PROD, BD_PBAR, BD_PBARL, BD_PBEAM, BD_PLOTEL
      USE SPRING_BUSH_MASS, ONLY      :  BD_CELAS1, BD_CELAS2, BD_CELAS3, BD_CELAS4, BD_PELAS, BD_CBUSH, BD_PBUSH, BD_CMASS1, BD_CMASS2, BD_CMASS3, BD_CMASS4, BD_PMASS, BD_CONM2
      USE SOLID_CARDS, ONLY           :  BD_CHEXA, BD_CPENTA, BD_CTETRA, BD_PSOLID
      USE GRID_COORDINATES, ONLY      :  BD_CORD, BD_GRID, BD_GRDSET, BD_SEQGP, BD_SPOINT, BD_SNORM
      USE SHELL_COMPOSITE_CARDS, ONLY :  BD_CQUAD, BD_CQUAD8, BD_CTRIA, BD_CSHEAR, BD_PSHEAR, BD_PSHEL, BD_PCOMP, BD_PCOMP1
      USE USER_ELEMENTS, ONLY         :  BD_CUSER1, BD_CUSERIN, BD_PUSER1, BD_PUSERIN
      USE DEBUG_CARDS, ONLY           :  BD_DEBUG
      USE EIGEN_NONLINEAR_CARDS, ONLY :  BD_EIGR, BD_EIGRL, BD_NLPARM
      USE BULK_DATA_LOADS, ONLY       :  BD_LOAD, BD_FORMOM, BD_GRAV, BD_PLOAD2, BD_PLOAD4, BD_RFORCE, BD_SLOAD
      USE MATERIAL_CARDS, ONLY        :  BD_MAT1, BD_MAT2, BD_MAT8, BD_MAT9
      USE CONSTRAINT_CARDS, ONLY      :  BD_SPC, BD_SPC1, BD_SPCADD, BD_MPC, BD_MPCADD
      USE PARAM_CARDS, ONLY           :  BD_PARAM
      USE PARTITION_VECTORS, ONLY     :  BD_PARVEC, BD_PARVEC1
      USE RIGID_ELEMENTS, ONLY        :  BD_RBAR, BD_RBE1, BD_RBE2, BD_RBE3, BD_RSPLINE
      USE BULK_DATA_TEMPERATURES, ONLY:  BD_TEMP, BD_TEMPD, BD_TEMPRP
      USE ALLOCATE_MODEL_STUF_Interface
      USE SORTING, ONLY               :  SORT_INT1

      END MODULE LOADB_USE_IFs
