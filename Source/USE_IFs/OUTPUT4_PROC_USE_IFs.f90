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

      MODULE OUTPUT4_PROC_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE OUTPUT4_PROC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_OU4_MAT_STATS
      USE OUTPUT4_PARTITIONING, ONLY  :  OU4_PARTVEC_PROC
      USE WRITE_PARTNd_MAT_HDRS_Interface
      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_FF, PARTITION_SS, PARTITION_SS_NTERM
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE WRITE_OU4_SPARSE_MAT_Interface
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
      USE WRITE_OU4_FULL_MAT_Interface
      USE OUTA_HERE_Interface

      END MODULE OUTPUT4_PROC_USE_IFs
