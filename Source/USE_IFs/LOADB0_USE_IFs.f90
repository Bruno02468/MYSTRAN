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

      MODULE LOADB0_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE LOADB0

      USE OURTIM_Interface
      USE OUTA_HERE_Interface
      USE FFIELD_Interface
      USE FFIELD2_Interface
      USE ROD_BAR_BEAM_CARDS, ONLY    :  BD_BAROR0, BD_BEAMOR0, BD_CBAR0
      USE SPRING_BUSH_MASS, ONLY      :  BD_CBUSH0
      USE SOLID_CARDS, ONLY           :  BD_CHEXA0, BD_CPENTA0, BD_CTETRA0
      USE SHELL_COMPOSITE_CARDS, ONLY :  BD_CQUAD0, BD_CQUAD80, BD_CTRIA0, BD_PCOMP0, BD_PCOMP10
      USE USER_ELEMENTS, ONLY         :  BD_CUSERIN0
      USE DEBUG_CARDS, ONLY           :  BD_DEBUG0
      USE GRID_COORDINATES, ONLY      :  BD_GRDSET0, BD_SPOINT0
      USE BULK_DATA_LOADS, ONLY       :  BD_LOAD0, BD_SLOAD0
      USE CONSTRAINT_CARDS, ONLY      :  BD_SPCADD0, BD_MPC0, BD_MPCADD0
      USE PARAM_CARDS, ONLY           :  BD_PARAM0
      USE RIGID_ELEMENTS, ONLY        :  BD_RBE30, BD_RSPLINE0

      END MODULE LOADB0_USE_IFs
