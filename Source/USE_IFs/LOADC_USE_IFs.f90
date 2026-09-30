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

      MODULE LOADC_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE LOADC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE OUTA_HERE_Interface
      USE REPLACE_TABS_W_BLANKS_Interface
      USE CSHIFT_Interface
      USE CASE_CONTROL_OUTPUTS, ONLY  :  CC_ACCE, CC_DISP, CC_ELDA, CC_ELFO, CC_ENFO, CC_GPFO, CC_MPCF, CC_OLOA, CC_SPCF, CC_STRE, CC_STRN
      USE CASE_CONTROL_METADATA, ONLY :  CC_ECHO, CC_LABE, CC_SUBT, CC_TITL
      USE CASE_CONTROL_SELECTORS, ONLY:  CC_LOAD, CC_METH, CC_MPC, CC_NLPARM, CC_SPC, CC_STATSUB, CC_TEMP
      USE CASE_CONTROL_SETS, ONLY     :  CC_SET, CC_SUBC

      END MODULE LOADC_USE_IFs
