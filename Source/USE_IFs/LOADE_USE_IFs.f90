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

      MODULE LOADE_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE LOADE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT, REPLACE_TABS_W_BLANKS
      USE EXECUTIVE_CONTROL_IO, ONLY  :  EC_DEBUG, EC_IN4FIL, EC_OUTPUT4, EC_PARTN
      USE BDF_SET_SYNTAX, ONLY        :  STOKEN
      USE BDF_FIELD_VALIDATION, ONLY  :  CRDERR, I4FLD

      END MODULE LOADE_USE_IFs
