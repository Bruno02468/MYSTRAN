! #################################################################################################################################
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
      SUBROUTINE CALC_K6ROT ()

! Builds the stiffness matrix for a spring connecting the drilling DOF to the translational DOFs of the adjacent nodes
! of the element. Adds that to the existing stiffness matrix KE in the element coordinate system.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  TYPE, ELGP, INTL_MID, XEL, SHELL_A, KE
      USE PARAMS, ONLY                :  K6ROT
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE SCONTR, ONLY                :  MAX_ORDER_GAUSS
      USE JACOBIAN, ONLY               :  JAC2D
      USE QUADRATURE, ONLY            :  ORDER_GAUSS
      USE VECTOR_GEOMETRY, ONLY       :  CROSS
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF_T

      IMPLICIT NONE

      REAL(DOUBLE)                    :: KROT(6*ELGP,6*ELGP)  ! Stifness matrix of K6ROT.
      REAL(DOUBLE)                    :: N(3)                 ! Normal at one grid point
      REAL(DOUBLE)                    :: X_PREV(3)            ! Coordinates of the previous grid point (n-1)
      REAL(DOUBLE)                    :: X_NEXT(3)            ! Coordinates of the next grid point (n+1)
      REAL(DOUBLE)                    :: TERM_PREV(3)
      REAL(DOUBLE)                    :: TERM_NEXT(3)
      REAL(DOUBLE)                    :: B(6*ELGP)            ! Strain-displacement matrix for one grid point's K6ROT spring "element"
      REAL(DOUBLE)                    :: STIFFNESS            ! Spring stiffness of one grid point's K6ROT "element"
      REAL(DOUBLE)                    :: AREA                 ! Elem area
      REAL(DOUBLE)                    :: DETJ                 ! An output from subr JAC2D4, called herein. Determinant of JAC
      INTEGER(LONG)                   :: I,J
      INTEGER(LONG)                   :: GP                   ! Element grid point number (1 to ELGP).
      INTEGER(LONG)                   :: GP_PREV
      INTEGER(LONG)                   :: GP_NEXT
      REAL(DOUBLE)                    :: JAC(2,2)             ! An output from subr JAC2D4, called herein. 2 x 2 Jacobian matrix.
      REAL(DOUBLE)                    :: JACI(2,2)            ! An output from subr JAC2D4, called herein. 2 x 2 Jacobian inverse.
      REAL(DOUBLE)                    :: X2E                  ! x coord of elem node 2
      REAL(DOUBLE)                    :: Y3E                  ! y coord of elem node 3
      REAL(DOUBLE)                    :: XSD(4)               ! Diffs in x coords of quad sides in local coords
      REAL(DOUBLE)                    :: YSD(4)               ! Diffs in y coords of quad sides in local coords
      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS) ! An output from subr ORDER, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: SSS(MAX_ORDER_GAUSS) ! An output from subr ORDER, called herein. Gauss abscissa's.


! **********************************************************************************************************************************


                                                     ! No K6ROT for shells that only use MID1.
      IF (INTL_MID(2) > 0) THEN

         AREA = ZERO

         IF ((TYPE(1:5) == "QUAD4")) THEN

            XSD(1) = XEL(1,1) - XEL(2,1)             ! x coord diffs (in local elem coords)
            XSD(2) = XEL(2,1) - XEL(3,1)
            XSD(3) = XEL(3,1) - XEL(4,1)
            XSD(4) = XEL(4,1) - XEL(1,1)

            YSD(1) = XEL(1,2) - XEL(2,2)             ! y coord diffs (in local elem coords)
            YSD(2) = XEL(2,2) - XEL(3,2)
            YSD(3) = XEL(3,2) - XEL(4,2)
            YSD(4) = XEL(4,2) - XEL(1,2)

            CALL ORDER_GAUSS ( 2, SSS, HHH )
            DO I=1,2
               DO J=1,2
                  CALL JAC2D ( SSS(I), SSS(J), XSD, YSD, 'N', JAC, JACI, DETJ )
                  AREA = AREA + HHH(I)*HHH(J)*DETJ
               ENDDO
            ENDDO

         ELSEIF (TYPE(1:5) == "TRIA3") THEN

            X2E  = XEL(2,1)
            Y3E  = XEL(3,2)
                                                     ! Actual area is half this but using this value
                                                     ! gives the same stiffness as MSC.
            AREA = X2E*Y3E

         ENDIF


         !According to:
         !https://www.dynalook.com/conferences/9th-european-ls-dyna-conference/drilling-rota
         !tion-constraint-for-shell-elements-in-implicit-and-explicit-analyses
         !but scaled to match MSC Nastran.

         ! SHELL_A(3,3) is membrane shear modulus times thickness
         STIFFNESS = 10.0**(-6.0) * K6ROT * SHELL_A(3,3) * ABS(AREA)

         ! Use the uniform element normal at all grid points.
         ! This might be supposed to be the shell normal but it hardly seems to make a difference and this way is simpler.
         N = [ZERO, ZERO, ONE]

         DO GP=1,ELGP

            B = ZERO

            ! The spring at one grid point is effectively a 3-node element with nodes and DOFs of:
            ! Node n:          tx, ty, tz, rx, ry, rz
            ! Node n+1 (next): tx, ty, tz
            ! Node n-1 (prev): tx, ty, tz

            GP_PREV = GP - 1
            IF (GP_PREV < 1) THEN
               GP_PREV = ELGP
            ENDIF

            GP_NEXT = GP + 1
            IF (GP_NEXT > ELGP) THEN
               GP_NEXT = 1
            ENDIF

            X_PREV = XEL(GP_PREV,1:3) - XEL(GP,1:3)
            X_NEXT = XEL(GP_NEXT,1:3) - XEL(GP,1:3)

            !Contribution of previous node's displacement
            !        - n × (x_n-1 - x_n)
            !ε_n  =  ------------------- * u_n-1
            !        2 * |x_n-1 - x_n|^2
            CALL CROSS(N, X_PREV, TERM_PREV)
            TERM_PREV = TERM_PREV * 1 / (2 * (X_PREV(1)**2 + X_PREV(2)**2 + X_PREV(3)**2) )

            B((GP_PREV - 1) * 6 + 1) = -TERM_PREV(1)
            B((GP_PREV - 1) * 6 + 2) = -TERM_PREV(2)
            B((GP_PREV - 1) * 6 + 3) = -TERM_PREV(3)

            !Contribution of next node's displacement
            !        - n × (x_n+1 - x_n)
            !ε_n += ------------------- * u_n+1
            !        2 * |x_n+1 - x_n|^2
            CALL CROSS(N, X_NEXT, TERM_NEXT)
            TERM_NEXT = TERM_NEXT * 1 / (2 * (X_NEXT(1)**2 + X_NEXT(2)**2 + X_NEXT(3)**2) )

            B((GP_NEXT - 1) * 6 + 1) = -TERM_NEXT(1)
            B((GP_NEXT - 1) * 6 + 2) = -TERM_NEXT(2)
            B((GP_NEXT - 1) * 6 + 3) = -TERM_NEXT(3)

            !Contribution of current node's displacement and rotation
            !          n × (x_n-1 - x_n)     n × (x_n+1 - x_n)
            !ε_n += ( ------------------- + ------------------- ) * u_n  +  n · r_n
            !         2 * |x_n-1 - x_n|^2   2 * |x_n+1 - x_n|^2
            B((GP - 1) * 6 + 1) = TERM_PREV(1) + TERM_NEXT(1)
            B((GP - 1) * 6 + 2) = TERM_PREV(2) + TERM_NEXT(2)
            B((GP - 1) * 6 + 3) = TERM_PREV(3) + TERM_NEXT(3)
            B((GP - 1) * 6 + 4) = N(1)
            B((GP - 1) * 6 + 5) = N(2)
            B((GP - 1) * 6 + 6) = N(3)

            ! stiffness * B' * B
            CALL MATMULT_FFF_T(B, B, 1, 6*ELGP, 6*ELGP, KROT)
            KROT = KROT * STIFFNESS

            KE(1:6*ELGP, 1:6*ELGP) = KE(1:6*ELGP, 1:6*ELGP) + KROT(1:6*ELGP, 1:6*ELGP)

         ENDDO

      ENDIF


! **********************************************************************************************************************************

      END SUBROUTINE CALC_K6ROT
