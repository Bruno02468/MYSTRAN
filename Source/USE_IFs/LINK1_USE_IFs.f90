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

      MODULE LINK1_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE LINK1

      USE TIME_INIT_Interface
      USE OURDAT_Interface
      USE OURTIM_Interface
      USE OUTA_HERE_Interface
      USE ALLOCATE_MODEL_STUF_Interface
      USE FILE_OPEN_Interface
      USE MPC_PROCESSING, ONLY        :  MPC_PROC
      USE FILE_CLOSE_Interface
      USE RIGID_ELEMENT_PROCESSING, ONLY:  RIGID_ELEM_PROC
      USE SPARSE_LOAD_CONSTRAINTS, ONLY:  SPARSE_PG, SPARSE_RMG
      USE FORCE_MOM_PROCESSING, ONLY  :  FORCE_MOM_PROC
      USE EPTL_Interface
      USE MASS_MATRIX_ASSEMBLY, ONLY  :  EMP0, EMP, MGGC_MASS_MATRIX, SPARSE_MGG
      USE ALLOCATE_EMS_ARRAYS_Interface
      USE ALLOCATE_L1_MGG_Interface
      USE DEALLOCATE_EMS_ARRAYS_Interface
      USE DEALLOCATE_L1_MGG_Interface
      USE DEALLOCATE_MODEL_STUF_Interface
      USE GRAV_PROCESSING, ONLY       :  GRAV_PROC
      USE RFORCE_PROCESSING, ONLY     :  RFORCE_PROC
      USE SLOAD_PROC_Interface
      USE STIFFNESS_MATRIX_ASSEMBLY, ONLY:  ESP0, ESP, SPARSE_KGG, SPARSE_KGGD
      USE ALLOCATE_STF_ARRAYS_Interface
      USE DEALLOCATE_IN4_FILES_Interface
      USE DEALLOCATE_STF_ARRAYS_Interface
      USE WRITE_DOF_TABLES_Interface
      USE ELEMENT_MODEL_INDEXING, ONLY:  ELSAVE
      USE CHK_ARRAY_ALLOC_STAT_Interface
      USE WRITE_ALLOC_MEM_TABLE_Interface
      USE WRITE_L1A_Interface
      USE FILE_INQUIRE_Interface

      END MODULE LINK1_USE_IFs
