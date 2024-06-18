! ##################################################################################################################################
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
  
! Generates shape functions for 2D elements AT THEIR MIDPOINTS (thus the M)
SUBROUTINE SHP2DM ( NMIDPT, CALLING_SUBR, IORD_MSG, IORZZZ, SSI, SSJ, WRT_BUG_THIS_TIME, PSH, DPSHG )

! The node numbering and axes convention are shown below with XI and ETA ranging from -1 to +1
!                         ETA
!                          |
!                          |
!                          |
!              . . . . . . 3 . . . . . .
!              .           |           .
!              .           |           .
!              .           |           .
!              4           ------------2----> XI
!              .                       .
!              .                       .
!              .                       .
!              . . . . . . 1 . . . . . .

   USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
   USE IOUNT1, ONLY                :  BUG, ERR, F04, F06, WRT_BUG, WRT_ERR, WRT_LOG
   USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_SHPJ_BIT, MEFE, FATAL_ERR
   USE TIMDAT, ONLY                :  TSEC
   USE SUBR_BEGEND_LEVELS, ONLY    :  SHP_BEGEND
   USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
   USE MODEL_STUF, ONLY            :  EID, EMG_IFE, ERR_SUB_NAM, NUM_EMG_FATAL_ERRS, TYPE

   USE SHP2DQ_USE_IFs

   IMPLICIT NONE

   CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SHP2DM'
   CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr (used for debug output)
   CHARACTER(LEN=*), INTENT(IN)    :: IORD_MSG          ! Character name of the integration order (used for debug output)
   CHARACTER(17*BYTE)              :: NAME(5)           ! Used for output annotation
   CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

   INTEGER(LONG), INTENT(IN)       :: NMIDPT            ! I index of Tying point (needed for some optional output)
   INTEGER(LONG), INTENT(IN)       :: IORZZZ            ! Integration order (used for debug output)
   INTEGER(LONG)                   :: I,J               ! DO loop indices
   INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = SHP_BEGEND

   REAL(DOUBLE) , INTENT(IN)       :: SSI               ! Tying point location component
   REAL(DOUBLE) , INTENT(IN)       :: SSJ               ! Tying point location component
   REAL(DOUBLE) , INTENT(OUT)      :: PSH(4)            ! Shape functions for all midpoints for this Tying point
   REAL(DOUBLE) , INTENT(OUT)      :: DPSHG(2,4)        ! Derivatives of PSH with respect to xi and eta.
   REAL(DOUBLE)                    :: XI(16), ET(16)    ! Elem node locations in isoparametric xi, eta coords
   REAL(DOUBLE)                    :: A1,A2,A3          ! Intermediate variables used in calculating outputs
   REAL(DOUBLE)                    :: B1,B2,B3          ! Intermediate variables used in calculating outputs
   REAL(DOUBLE)                    :: XI2               ! Squares of xi coords
   REAL(DOUBLE)                    :: ET2               ! Squares of eta coords

   ! **********************************************************************************************************************************
   IF (WRT_LOG >= SUBR_BEGEND) THEN
      CALL OURTIM
      WRITE(F04,9001) SUBR_NAME,TSEC, WRT_BUG_THIS_TIME, WRT_BUG(7), WRT_BUG(8), WRT_BUG(9)
   9001 FORMAT(1X,A,' BEGN ',F10.3, 3X, A1, 3(I3))
   ENDIF
   
   ! **********************************************************************************************************************************
   ! Initialize outputs

   DO I=1,4
      PSH(I) = ZERO
   ENDDO

   DO I=1,2
      DO J=1,4
         DPSHG(I,J) = ZERO
      ENDDO
   ENDDO

   ! grid points
   XI(1) = -ONE
   XI(2) =  ONE
   XI(3) =  ONE
   XI(4) = -ONE
   
   ! midpoints
   !XI(1) = ZERO
   !XI(2) = ONE
   !XI(3) = ZERO
   !XI(4) = -ONE

   ! grid points
   ET(1) = -ONE
   ET(2) = -ONE
   ET(3) =  ONE
   ET(4) =  ONE

   ! midpoints
   !ET(1) = -ONE
   !ET(2) = ZERO
   !ET(3) = ONE
   !ET(4) = ZERO


   DO I=1,4
      PSH(I) = (ONE + SSI*XI(I))*(ONE + SSJ*ET(I))/FOUR
      DPSHG(1,I) = XI(I)*(ONE + SSJ*ET(I))/FOUR
      DPSHG(2,I) = ET(I)*(ONE + SSI*XI(I))/FOUR
   ENDDO 

END SUBROUTINE SHP2DM
