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
 
      SUBROUTINE NEXTC ( CARD, ICONTINUE, IERR )
 
      ! Looks for a Bulk Data continuation card belonging to a parent card.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_LOG, ERR, F04, F06, IN1, INFILE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, ECHO, FATAL_ERR, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC
      USE SUBR_BEGEND_LEVELS, ONLY    :  NEXTC_BEGEND
 
      USE NEXTC_USE_IFs

      IMPLICIT NONE
 
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'NEXTC'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A MYSTRAN data card
      CHARACTER(LEN=LEN(CARD))        :: CARD_IN           ! Version of CARD read here
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10), JCARD0(10) ! 10 fields of 8 characters of CARD
      CHARACTER(24*BYTE)              :: MESSAG            ! Message for output error purposes
      CHARACTER(LEN(JCARD))           :: NEWTAG            ! Field 1  of cont   card
      CHARACTER(LEN(JCARD))           :: OLDTAG            ! Field 10 of parent card
      CHARACTER(LEN=LEN(CARD))        :: TCARD             ! Temporary version of CARD
 
      INTEGER(LONG), INTENT(OUT)      :: ICONTINUE         ! =1 if next card is current card's continuation or =0 if not
      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator from subr FFIELD, called herein
      INTEGER(LONG)                   :: COMMENT_COL       ! Col on CARD where a comment begins (if one exists)
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from READ
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR  
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file. Input to subr READERR
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = NEXTC_BEGEND
 
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
      ! Initialize error indicator and ICONTINUE
      IERR = 0
      ICONTINUE = 0

      ! Make units for writing errors the error file and output file
      OUNT(1) = ERR
      OUNT(2) = F06

      ! Make JCARD for parent CARD
      ! split the line (CARD) into fields (JCARD)
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      
      ! copy jcard to jcard0 in case we have an error
      DO I=1,10
          JCARD0(I) = JCARD(I)
      ENDDO

      !------------------------
      ! Read next card.

      ! OLDTAG is field 10 of the current card coming into this subr
      OLDTAG = JCARD(10)
      MESSAG = 'BULK DATA CARD          '

      ! Read next card
      CALL READ_BDF_LINE(IN1, IOCHK, TCARD)
      CARD_IN = TCARD
      IF (IOCHK /= 0) THEN
         REC_NO = -99
         CALL READERR ( IOCHK, INFILE, MESSAG, REC_NO, OUNT, 'Y' )
         FATAL_ERR = FATAL_ERR + 1
      ENDIF
      
      ! we (maybe) found a comment line inside the card
      ! loop until we find a non-comment line using the same
      ! block as above
      !
      ! NOTE: comment lines must start with a $, so a space is
      !       not considered a comment
      DO WHILE (TCARD(1:1) == '$') 
         ! Read next card line
         CALL READ_BDF_LINE(IN1, IOCHK, TCARD)
         CARD_IN = TCARD
         IF (IOCHK /= 0) THEN
            REC_NO = -99
            CALL READERR ( IOCHK, INFILE, MESSAG, REC_NO, OUNT, 'Y' )
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO

      ! Remove any comments within the CARD by deleting everything 
      ! from $ on (after col 1)
      COMMENT_COL = 1
      DO I=2,BD_ENTRY_LEN
         IF (TCARD(I:I) == '$') THEN
            COMMENT_COL = I
            EXIT
         ENDIF
      ENDDO

      IF (COMMENT_COL > 1) THEN
         TCARD(COMMENT_COL:) = ' '
      ENDIF

      ! Make JCARD for TCARD above and get FFIELD to left adjust and
      ! fix-field it (if necessary).
      !
      ! get a flag (ICONTINUE) that defines if we need to read the next line
      ! we do because we have OLDTAG

      ! we know know that the line doesn't start with a $, but it
      ! can be a different card (OLDTAG /= NEWTAG)
      !
      CALL FFIELD ( TCARD, IERR )

      ! split the continuation line (TCARD) into fields (JCARD)
      CALL MKJCARD ( SUBR_NAME, TCARD, JCARD )
      NEWTAG = JCARD(1)
      
      !WRITE(ERR,*) 'tags',OLDTAG,NEWTAG
      
      IF (NEWTAG == OLDTAG) THEN
         ICONTINUE = 1
      ELSE IF ((OLDTAG(1:1) == '+') .AND. (NEWTAG(1:1) == ' ') .AND. (OLDTAG(2:8) == NEWTAG(2:8))) THEN
         ICONTINUE = 1
      ELSE IF ((OLDTAG(1:1) == ' ') .AND. (NEWTAG(1:1) == '+') .AND. (OLDTAG(2:8) == NEWTAG(2:8))) THEN
         ICONTINUE = 1
      ELSE IF ((NEWTAG(1:1) /= ' ') .AND. (NEWTAG(1:1) /= '+') .AND. (NEWTAG(1:1) /= '$')) THEN
         ! different card type (e.g., LOAD -> FORCE
         BACKSPACE(IN1)
      ELSE
         ! can\t find the continuation marker.  FATAL :)
         BACKSPACE(IN1)
         WRITE(F06,102) OLDTAG
         WRITE(ERR,102) OLDTAG
         WRITE(F06,103)
         WRITE(ERR,103)
         WRITE(F06,104) 'FIELDS1:', JCARD0
         WRITE(ERR,104) 'FIELDS1:', JCARD0
         WRITE(F06,104) 'FIELDS2:', JCARD
         WRITE(ERR,104) 'FIELDS2:', JCARD
         FLUSH(F06)
         FLUSH(ERR)
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE('Y')  ! FATAL error
         RETURN
      ENDIF
      CARD = TCARD
      IF (ECHO(1:4) /= 'NONE') THEN
          WRITE(F06,'(A)') CARD_IN
      ENDIF

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

  ! missing continuation fatal
  102 FORMAT(' *FATAL: COULD NOT FIND A CONTINUATION MARKER FOR ', A)
  103 FORMAT('         CHECK THAT THE CONTINUATION STARTS AT THE BEGINNING OF FIELD 10')

  ! missing continuation; name and 10 fields
  104 FORMAT('         ', A, ' 1: ', A, ' 2:', A, ' 3:', A, ' 4:', A, &
             ' 5:', A, ' 6:', A, ' 7:', A, ' 8:', A, ' 9:', A, ' 10:', A)
! **********************************************************************************************************************************

      END SUBROUTINE NEXTC
