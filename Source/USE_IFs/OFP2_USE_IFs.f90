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

      MODULE OFP2_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE OFP2

      USE OURTIM_Interface
      USE ALLOCATE_COL_VEC_Interface
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE MATMULT_SFF_Interface
      USE DEALLOCATE_COL_VEC_Interface
      USE GET_ARRAY_ROW_NUM_Interface
      USE GET_GRID_NUM_COMPS_Interface
      USE CHK_OGEL_ZEROS_Interface
      USE GRID_OUTPUT_WRITERS, ONLY   :  WRITE_FEMAP_GRID_VECS, WRITE_GRD_PRT_OUTPUTS
      USE CONVERT_VEC_COORD_SYS_Interface
      USE MERGE_COL_VECS_Interface
      USE OUTA_HERE_Interface

      END MODULE OFP2_USE_IFs
