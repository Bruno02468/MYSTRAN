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

      MODULE MYSTRAN_USE_IFs

! USE Interface statements for all subroutines called by PROGRAM MYSTRAN

      USE DATE_TIME_UTILS, ONLY       :  OURDAT, OURTIM, TIME_INIT
      USE INI_FILE, ONLY            :  READ_INI
      USE READ_INPUT_FILE_NAME_Interface
      USE FILE_LIFECYCLE, ONLY        :  CLOSE_LIJFILES, CLOSE_OUTFILES, FILE_CLOSE, FILE_INQUIRE, FILE_OPEN, WRITE_FILNAM
      USE IS_THIS_A_RESTART_Interface
      USE MYSTRAN_FILES_Interface
      USE PROCESS_INCLUDE_FILES_Interface
      USE LOADE0_Interface
      USE TEMP_FILE_READERS, ONLY     :  READ_L1A
      USE LINK0_LINK1_MOD, ONLY       :  LINK0, LINK1
      USE LINK2_MOD, ONLY             :  LINK2
      USE LINK3_MOD, ONLY             :  LINK3
      USE LINK4_MOD, ONLY             :  LINK4
      USE LINK6_MOD, ONLY             :  LINK6
      USE RIGID_BODY_STORAGE_LIFECYCLE, ONLY:  DEALLOCATE_RBGLOBAL
      USE LINK5_MOD, ONLY             :  LINK5
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE LINK9_MOD, ONLY             :  LINK9
      USE RESTART_FILE_IO, ONLY       :  RESTART_DATA_FOR_L3
      USE NONLINEAR_PARAMETER_LIFECYCLE, ONLY:  DEALLOCATE_NL_PARAMS
      USE VECTOR_METRICS, ONLY        :  VECTOR_NORM
      USE PRINT_BUILD_INFO_Interface
      USE READ_CL_Interface
      USE SET_BLAS_THREADS_Interface

      END MODULE MYSTRAN_USE_IFs
