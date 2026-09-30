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

      MODULE SEQ_PROC_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE SEQ_PROC

      USE OURTIM_Interface
      USE GET_ARRAY_ROW_NUM_Interface
      USE OUTA_HERE_Interface
      USE SORTING, ONLY               :  SORT_INT1, SORT_INT1_REAL1, SORT_INT2, SORT_INT2_REAL1, SORT_REAL1_INT1
      USE OPNERR_Interface
      USE READERR_Interface
      USE STMERR_Interface
      USE DEALLOCATE_MODEL_STUF_Interface
      USE ALLOCATE_MODEL_STUF_Interface
      USE BDF_CARD_CONTINUATIONS, ONLY:  MKCARD, MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  LEFT_ADJ_BDFLD
      USE GRID_COORDINATES, ONLY      :  BD_SEQGP
      USE WRITE_INTEGER_VEC_Interface
      USE FILE_CLOSE_Interface

      END MODULE SEQ_PROC_USE_IFs
