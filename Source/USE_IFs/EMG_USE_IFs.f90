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

      MODULE EMG_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE EMG

      USE COMPOSITE_SHELL_PREPARATION, ONLY:  IS_ELEM_PCOMP_PROPS, SHELL_ABD_MATRICES
      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE ELEMENT_DATA_GATHERING, ONLY:  ELMDAT1, ELMDAT2
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE ELEMENT_GEOMETRY_PREPARATION, ONLY:  ELMGM1, ELMGM2, ELMGM3
      USE ELEMENT_LOOKUPS, ONLY       :  GET_MATANGLE_FROM_CID
      USE MATERIAL_PROPERTIES, ONLY   :  MATERIAL_PROPS_2D, MATERIAL_PROPS_3D
      USE MATERIAL_TRANSFORMATIONS, ONLY:  ROT_AXES_MATL_TO_LOC
      USE ELEMENT_DIAGNOSTICS, ONLY   :  ELMOFF, ELMOUT
      USE SPRING_ELEMENTS, ONLY       :  BUSH, ELAS1
      USE LINE_ELEMENTS, ONLY         :  BREL1
      USE TRIANGULAR_SHELL_ELEMENTS, ONLY:  TREL1
      USE QUADRILATERAL_ELEMENT_DISPATCH, ONLY:  QDEL1
      USE SOLID_ELEMENTS, ONLY        :  HEXA, PENTA, TETRA
      USE USER_DEFINED_ELEMENTS, ONLY :  KUSER1, USERIN
      USE MITC8_MOD, ONLY             :  MITC8

      END MODULE EMG_USE_IFs
