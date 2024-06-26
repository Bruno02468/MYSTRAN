! ###############################################################################################################################
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

MODULE BSMITC_Interface

   INTERFACE

      SUBROUTINE BSMITC ( PSH, DPSHX, DNXSHX, DNYSHX, NMIDPT, MESSAG, WRT_BUG_THIS_TIME, BS )

   
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, F04, WRT_BUG, WRT_LOG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_BMAT_BIT, ELDT_BUG_BCHK_BIT
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE SUBR_BEGEND_LEVELS, ONLY    :  BSMITC_BEGEND
      USE MODEL_STUF, ONLY            :  EID, TYPE, XEB, XEL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
   
      IMPLICIT NONE
   
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Messag to print out if BCHECK is run
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG), INTENT(IN)       :: NMIDPT             ! I index of tying point (needed for some optional output)
      INTEGER(LONG), PARAMETER        :: NR        = 2     ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: NC        = 12    ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = BSMITC_BEGEND
   
      REAL(DOUBLE) , INTENT(IN)       :: PSH(4)            ! 4 node bilinear isopar interp functions (used for bending)
      REAL(DOUBLE) , INTENT(IN)       :: DPSHX(2,4)        ! Derivatives of PSH shape functions wrt x and y
      REAL(DOUBLE) , INTENT(IN)       :: DNXSHX(2,4)       ! Derivatives of constrained interpolations NX wrt x, y
      REAL(DOUBLE) , INTENT(IN)       :: DNYSHX(2,4)       ! Derivatives of constrained interpolations NY wrt x, y
      REAL(DOUBLE) , INTENT(OUT)      :: BS(2,12)          ! Output strain-displ matrix for this elem
      END SUBROUTINE BSMITC

   END INTERFACE
   
END MODULE BSMITC_Interface
   
   