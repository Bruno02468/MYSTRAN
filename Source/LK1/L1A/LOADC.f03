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

      SUBROUTINE LOADC

      ! LOADC reads in the CASE CONTROL DECK
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUGOUT, ERR, F06, IN1, WRT_ERR
      USE SCONTR, ONLY                :  CC_ENTRY_LEN, ENFORCED, FATAL_ERR, WARN_ERR, NSUB, NTSUB,                                 &
                                         NUM_BUCKLING_SUBS, PROG_NAME, RESTART, SOL_NAME, CC_CMD_DESCRIBERS
      USE PARAMS, ONLY                :  SUPINFO, SUPWARN
      USE MODEL_STUF, ONLY            :  CC_EIGR_SID, CC_EIGR_SID_SUB, CC_EIGR_SID_DECK, CC_STATSUB_DECK, CC_STATSUB_SUB,          &
                                         IS_BUCKLING_SUBCASE, IS_MODES_SUBCASE,                                                    &
                                         MEFFMASS_CALC, MPCSET, MPCSETS, MPFACTOR_CALC, SCNUM, SPCSET, SPCSETS, SUBLOD,            &
                                         SC_STRE, SC_STRN, SC_ELFE, SC_ELFN
      USE MODEL_STUF, ONLY            :  EIG_PARAMS,                                                                               &
                                         EIG_COMP, EIG_CRIT, EIG_FRQ1, EIG_FRQ2, EIG_GRID, EIG_LANCZOS_NEV_DELT, EIG_METH,         &
                                         EIG_MSGLVL, EIG_LAP_MAT_TYPE, EIG_MODE, EIG_N1, EIG_N2, EIG_NCVFACL, EIG_NORM, EIG_SID,   &
                                         EIG_SIGMA, EIG_VECS
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_LOC, STRE_LOC, FORC_LOC

      USE LOADC_USE_IFs
      USE TEXT_FIELD_UTILS, ONLY      :  TO_UPPER
      
      IMPLICIT NONE

      CHARACTER( 1*BYTE)              :: DOLLAR_WARN       ! Indicator of whether there was a $ sign in col 1
      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD              ! Case Control card
      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD1             ! CARD shifted to begin in col 1
      CHARACTER(12*BYTE)              :: DECK_NAME   = 'CASE CONTROL'
      CHARACTER(10*BYTE), PARAMETER   :: END_CARD    = 'BEGIN BULK'

      INTEGER(LONG)                   :: CHAR_COL          ! Column number on CARD where character CHAR is found
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! Error indicator. If CHAR not found, IERR set to 1
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a Case Control card from unit IN1

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: QUAD4_LOC


! **********************************************************************************************************************************

      ! Process CASE CONTROL DECK

outer:DO
         DOLLAR_WARN = 'N'
         CALL READ_BDF_LINE(IN1, IOCHK, CARD)

         ! Quit if EOF/EOR occurs during read
         IF (IOCHK < 0) THEN
            WRITE(ERR,1011) END_CARD
            WRITE(F06,1011) END_CARD
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         ! Check if error occurs during read.
         IF (IOCHK > 0) THEN
            WRITE(ERR,1010) DECK_NAME
            WRITE(F06,1010) DECK_NAME
            WRITE(F06,'(A)') CARD
            FATAL_ERR = FATAL_ERR + 1
            CYCLE outer
         ENDIF

         WRITE(F06,101) CARD

         ! Replace all tab characters with a white space
         CALL REPLACE_TABS_W_BLANKS ( CARD )

         ! Shift card so that it begins in col 1
         CALL CSHIFT ( CARD, ' ', CARD1, CHAR_COL, IERR )

         ! Check for CASE CONTROL cards. Exit loop on 'BEGIN BULK'
         IF      (CARD1(1:4) == 'ACCE'    ) THEN
            CALL CC_ACCE ( CARD1 )

         ELSE IF(CARD1(1:10) == 'BEGIN BULK') THEN
            IF (DOLLAR_WARN == 'Y') THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1199) CARD1
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,1199) CARD1
               ENDIF
            ENDIF
            EXIT outer

         ELSE IF((CARD1(1:4) == 'DISP'    ) .OR.  (CARD1(1:6) == 'VECTOR'  )) THEN
            CALL CC_DISP   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'ECHO'    ) THEN
            CALL CC_ECHO   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'ELDA'    ) THEN
            CALL CC_ELDA   ( CARD1 )
            BUGOUT = 'Y'

         ELSE IF((CARD1(1:7) == 'ELFORCE') .OR. (CARD1(1:4) == 'ELFO').OR. (CARD1(1:5) == 'FORCE')) THEN
            CALL CC_ELFO   ( CARD1 )

         ELSE IF (CARD1(1:8) == 'ENFORCED') THEN
            CALL CC_ENFO   ( CARD1 )
            ENFORCED = 'Y'

         ELSE IF (CARD1(1:4) == 'GPFO'    ) THEN
            CALL CC_GPFO   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'LABE'    ) THEN
            CALL CC_LABE   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'LOAD'    ) THEN
            CALL CC_LOAD   ( CARD1 )

         ELSE IF (CARD1(1:8) == 'MEFFMASS') THEN
            IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
               MEFFMASS_CALC = 'Y'
            ENDIF

         ELSE IF (CARD1(1:4) == 'METH'    ) THEN
            CALL CC_METH   ( CARD1 )

         ELSE IF((CARD1(1:3) == 'MPC'     ) .AND. (CARD1(1:4) /= 'MPCF'    )) THEN
            CALL CC_MPC    ( CARD1 )

         ELSE IF (CARD1(1:4) == 'MPCF'    ) THEN
            CALL CC_MPCF   ( CARD1 )

         ELSE IF (CARD1(1:8) == 'MPFACTOR') THEN
            IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
               MPFACTOR_CALC = 'Y'
            ENDIF

         ELSE IF (CARD1(1:6) == 'NLPARM'  ) THEN
            CALL CC_NLPARM ( CARD1 )

         ELSE IF (CARD1(1:4) == 'OLOA'    ) THEN
            CALL CC_OLOA   ( CARD1 )

         ELSE IF (CARD1(1:6) == 'OUTPUT' ) THEN            ! Normal OUTPUT entry is OK. Ones like OUTPUT(PLOT) we end CC processing.
            IF (INDEX(CARD,"(") > 0) THEN                  ! If we find "(" in an OUTPUT CC entry it indicates, e.g., OUTPUT(PLOT),
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,902) CARD
               WRITE(F06,902) CARD

               ! so read all entries up to BEGIN BULK and then exit loop for CC entries
inner:         DO
                  READ(IN1,101) CARD
                  CARD = TO_UPPER(CARD)
                  IF (CARD(1:10) == 'BEGIN BULK') THEN
                     WRITE(F06,101) CARD
                     EXIT outer
                  ELSE
                     CYCLE inner
                  ENDIF
               ENDDO inner
            ELSE                                           ! The OUTPUT entry was a normal one (with no "(PLOT)", etc delimiter,
               CYCLE outer                                 ! so continue processing CC entries
            ENDIF

         ELSE IF (CARD1(1:3) == 'SET'     ) THEN
            CALL CC_SET    ( CARD1 )

         ELSE IF((CARD1(1:3) == 'SPC'     ) .AND. (CARD1(1:4) /= 'SPCF'    )) THEN
            CALL CC_SPC    ( CARD1 )

         ELSE IF (CARD1(1:4) == 'SPCF'    ) THEN
            CALL CC_SPCF   ( CARD1 )

         ELSE IF (CARD1(1:7) == 'STATSUB') THEN
            CALL CC_STATSUB ( CARD1 )

         ELSE IF((CARD1(1:4) == 'STRA'    ) .OR.  (CARD1(1:4) == 'STRN'    ) .OR.  (CARD1(1:8) == 'ELSTRAIN')) THEN
            CALL CC_STRN   ( CARD1 )

         ELSE IF((CARD1(1:4) == 'STRE'    ) .OR.  (CARD1(1:4) == 'STRS'    ) .OR.  (CARD1(1:8) == 'ELSTRESS')) THEN
            CALL CC_STRE   ( CARD1 )

         ELSE IF (CARD1(1:8) == 'SUBCASE ') THEN
            CALL CC_SUBC   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'SUBT'    ) THEN
            CALL CC_SUBT   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'TEMP'    ) THEN
            CALL CC_TEMP   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'TITL'    ) THEN
            CALL CC_TITL   ( CARD1 )

         ELSE IF (CARD(1:1) == '$') THEN
            DO I=IACHAR('A'),IACHAR('Z')
               IF (CARD(2:2) == ACHAR(I)) THEN
                  DOLLAR_WARN = 'Y'
               ENDIF
            ENDDO

         ELSE IF (CARD(1:CC_ENTRY_LEN) == ' ') THEN
            CYCLE outer

         ELSE                                              ! Card is not recognized.
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) CARD
            WRITE(ERR,9993) PROG_NAME
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9993) PROG_NAME
            ENDIF

         ENDIF

      ENDDO outer

                                                           ! Assign the same CENTER/CORNER location to all of FORCE,
                                                           ! STRESS, and STRAIN outputs according to the priority rules.
                                                           ! Because NONE is stored as 0 in SC_STRE, etc., this treats 
                                                           ! STRESS(CORNER) = NONE as if STRESS wasn't defined at
                                                           ! all and defaults to CENTER.
      IF       (SC_STRE(1) /= 0) THEN
         QUAD4_LOC = STRE_LOC
      ELSE IF  (SC_STRN(1) /= 0) THEN
         QUAD4_LOC = STRN_LOC
      ELSE IF  (SC_ELFE(1) /= 0 .OR. SC_ELFN(1) /= 0) THEN
         QUAD4_LOC = FORC_LOC
      ELSE
         QUAD4_LOC = 'CENTER'
      ENDIF
      STRE_LOC = QUAD4_LOC
      STRN_LOC = QUAD4_LOC
      FORC_LOC = QUAD4_LOC
                                                           ! From here on, STRE_LOC = STRN_LOC = FORC_LOC.

      IF (STRN_LOC /= 'CENTER  ') THEN
         WRITE(ERR,1016) STRN_LOC, 'STRAIN'
         IF (SUPINFO == 'N') THEN
            WRITE(F06,1016) STRN_LOC, 'STRAIN'
         ENDIF
      ENDIF
      IF (STRE_LOC /= 'CENTER  ') THEN
         WRITE(ERR,1016) STRE_LOC, 'STRESS'
         IF (SUPINFO == 'N') THEN
            WRITE(F06,1016) STRE_LOC, 'STRESS'
         ENDIF
      ENDIF
      IF (FORC_LOC /= 'CENTER  ') THEN
         WRITE(ERR,1016) FORC_LOC, 'FORCE'
         IF (SUPINFO == 'N') THEN
            WRITE(F06,1016) FORC_LOC, 'FORCE'
         ENDIF
      ENDIF

      IF (NSUB == 0) THEN                                  ! There was no SUBCASE entry in Case Control so set = 1
         NSUB     = 1
         SCNUM(1) = 1
      ENDIF

      ! Propagate the deck-level METHOD default down to every subcase that does not have its own METHOD.
      ! This must happen *after* the NSUB==0 -> NSUB=1 fixup above, so that an "implicit" single subcase still picks
      ! up the METHOD declared at the deck level. The deck default is also recorded as the SID for any subcase that
      ! omitted METHOD. For SOL types that do not consume eigenvalue extraction the arrays simply stay zero.
      IF (ALLOCATED(CC_EIGR_SID_SUB)) THEN
         DO I=1,NSUB
            IF ((CC_EIGR_SID_SUB(I) == 0) .AND. (CC_EIGR_SID_DECK /= 0)) THEN
               CC_EIGR_SID_SUB(I)  = CC_EIGR_SID_DECK
               IS_MODES_SUBCASE(I) = 'Y'
            ENDIF
         ENDDO
         ! Keep the legacy scalar in sync: prefer the deck default; otherwise the first non-zero per-subcase SID.
         IF (CC_EIGR_SID_DECK /= 0) THEN
            CC_EIGR_SID = CC_EIGR_SID_DECK
         ELSE
            DO I=1,NSUB
               IF (CC_EIGR_SID_SUB(I) /= 0) THEN
                  CC_EIGR_SID = CC_EIGR_SID_SUB(I)
                  EXIT
               ENDIF
            ENDDO
         ENDIF

         ! Final fallback: for any subcase whose EIG_PARAMS slot is still empty (could happen when NSUB was 0 at the
         ! time BD_EIGR/BD_EIGRL ran -- e.g. a deck with a deck-level METHOD and no SUBCASE cards), copy the values
         ! from the legacy EIG_* scalars. After BD_EIGR/BD_EIGRL has WRITTEN L1M for the canonical (scalar-match) card,
         ! those scalars hold that card's data, which is exactly what an inherited subcase should use.
         IF (ALLOCATED(EIG_PARAMS)) THEN
            DO I=1,NSUB
               IF ((CC_EIGR_SID_SUB(I) /= 0) .AND. (EIG_PARAMS(I)%SID == 0)) THEN
                  EIG_PARAMS(I)%METHOD            = EIG_METH
                  EIG_PARAMS(I)%NORM              = EIG_NORM
                  EIG_PARAMS(I)%LAP_MAT_TYPE      = EIG_LAP_MAT_TYPE
                  EIG_PARAMS(I)%VECS              = EIG_VECS
                  EIG_PARAMS(I)%SID               = CC_EIGR_SID_SUB(I)
                  EIG_PARAMS(I)%N1                = EIG_N1
                  EIG_PARAMS(I)%N2                = EIG_N2
                  EIG_PARAMS(I)%COMP              = EIG_COMP
                  EIG_PARAMS(I)%GRID              = EIG_GRID
                  EIG_PARAMS(I)%LANCZOS_NEV_DELT  = EIG_LANCZOS_NEV_DELT
                  EIG_PARAMS(I)%MODE              = EIG_MODE
                  EIG_PARAMS(I)%MSGLVL            = EIG_MSGLVL
                  EIG_PARAMS(I)%NCVFACL           = EIG_NCVFACL
                  EIG_PARAMS(I)%CRIT              = EIG_CRIT
                  EIG_PARAMS(I)%FRQ1              = EIG_FRQ1
                  EIG_PARAMS(I)%FRQ2              = EIG_FRQ2
                  EIG_PARAMS(I)%SIGMA             = EIG_SIGMA
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      ! If SOL is modes or CB or buckilng, then a METH card should have been found in Case Control
      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN_CB_MODEL') .OR. (SOL_NAME(1:8) == 'BUCKLING')) THEN
         IF (CC_EIGR_SID == 0) THEN
             WRITE(ERR,1004)
             WRITE(F06,1004)
             FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

     ! For SOL 105 buckling, resolve STATSUB(PRELOAD) references for each buckling subcase. A buckling subcase is one
     ! that has a resolved METHOD (and therefore IS_MODES_SUBCASE(I)=='Y'). The remaining subcases supply linear-static
     ! preloads. STATSUB(PRELOAD)=n names the external SUBCASE id of the static subcase whose displacement field is
     ! used to assemble KGGD for the buckling eigenproblem. If no STATSUB is given anywhere, fall back to the legacy
     ! behavior: the first non-buckling subcase that carries a LOAD or TEMP request acts as the preload source.
      IF (SOL_NAME == 'BUCKLING') THEN

         ! Tag every modes-subcase as a buckling-subcase, and count them.
         NUM_BUCKLING_SUBS = 0
         IF (ALLOCATED(IS_BUCKLING_SUBCASE) .AND. ALLOCATED(IS_MODES_SUBCASE)) THEN
            DO I = 1, NSUB
               IF (IS_MODES_SUBCASE(I) == 'Y') THEN
                  IS_BUCKLING_SUBCASE(I) = 'Y'
                  NUM_BUCKLING_SUBS = NUM_BUCKLING_SUBS + 1
               ENDIF
            ENDDO
         ENDIF

         IF (NUM_BUCKLING_SUBS == 0) THEN                  ! no buckling subcase -> nothing to solve
            WRITE(ERR,1101)
            WRITE(F06,1101)
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         ! Confirm at least one non-buckling subcase carries a load (mechanical or thermal). Without one there is no
         ! linear-static preload to drive KGGD.
         BLOCK
            LOGICAL :: HAS_STATIC_LOAD
            INTEGER(LONG) :: STATSUB_REQ, RESOLVED_IDX, K
            HAS_STATIC_LOAD = .FALSE.
            DO I = 1, NSUB
               IF ((IS_BUCKLING_SUBCASE(I) == 'N') .AND. ((SUBLOD(I,1) /= 0) .OR. (SUBLOD(I,2) /= 0))) THEN
                  HAS_STATIC_LOAD = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (.NOT. HAS_STATIC_LOAD) THEN
               WRITE(ERR,1101)
               WRITE(F06,1101)
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

            ! Resolve STATSUB(PRELOAD) for each buckling subcase. Per-subcase value beats deck-level which beats the
            ! legacy fallback (first static subcase that carries a load).
            DO I = 1, NSUB
               IF (IS_BUCKLING_SUBCASE(I) /= 'Y') CYCLE
               STATSUB_REQ = 0
               IF (ALLOCATED(CC_STATSUB_SUB)) STATSUB_REQ = CC_STATSUB_SUB(I)
               IF (STATSUB_REQ == 0) STATSUB_REQ = CC_STATSUB_DECK

               RESOLVED_IDX = 0
               IF (STATSUB_REQ == 0) THEN
                  ! Legacy fallback: first static subcase with LOAD or TEMP.
                  DO K = 1, NSUB
                     IF ((IS_BUCKLING_SUBCASE(K) == 'N') .AND. ((SUBLOD(K,1) /= 0) .OR. (SUBLOD(K,2) /= 0))) THEN
                        RESOLVED_IDX = K
                        EXIT
                     ENDIF
                  ENDDO
                  IF (RESOLVED_IDX == 0) THEN
                     WRITE(ERR,1103) I
                     WRITE(F06,1103) I
                     FATAL_ERR = FATAL_ERR + 1
                  ENDIF
               ELSE
                  ! Look up external subcase id (STATSUB_REQ) in SCNUM(:) to find its internal index.
                  DO K = 1, NSUB
                     IF (SCNUM(K) == STATSUB_REQ) THEN
                        RESOLVED_IDX = K
                        EXIT
                     ENDIF
                  ENDDO
                  IF (RESOLVED_IDX == 0) THEN
                     WRITE(ERR,1104) STATSUB_REQ, I
                     WRITE(F06,1104) STATSUB_REQ, I
                     FATAL_ERR = FATAL_ERR + 1
                  ELSE IF (IS_BUCKLING_SUBCASE(RESOLVED_IDX) == 'Y') THEN
                     WRITE(ERR,1105) STATSUB_REQ, I
                     WRITE(F06,1105) STATSUB_REQ, I
                     FATAL_ERR = FATAL_ERR + 1
                     RESOLVED_IDX = 0
                  ELSE IF ((SUBLOD(RESOLVED_IDX,1) == 0) .AND. (SUBLOD(RESOLVED_IDX,2) == 0)) THEN
                     WRITE(ERR,1106) STATSUB_REQ, I
                     WRITE(F06,1106) STATSUB_REQ, I
                     FATAL_ERR = FATAL_ERR + 1
                  ENDIF
               ENDIF

               IF (ALLOCATED(EIG_PARAMS) .AND. (RESOLVED_IDX > 0)) THEN
                  EIG_PARAMS(I)%STATSUB_REF = RESOLVED_IDX
               ENDIF
            ENDDO
         END BLOCK

      ENDIF

      ! Make sure that NTSUB <= NSUB (if user had more TEMP cards in C.C.
      ! than subcases, this will give problems with size of TCASE2)
      IF (NTSUB > NSUB) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1803) NTSUB, NSUB
         WRITE(F06,1803) NTSUB, NSUB
      ENDIF


     ! Make sure that, if there are more than 1 of SPC or MPC set requests
     ! in Case Control, that all (nonzero) are the same. Otherwise issue error.
      SPCSET = SPCSETS(1)
      DO I=2,NSUB
         IF (SPCSETS(I) /= 0) THEN
            IF (SPCSETS(I) /= SPCSET) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1830) 'SPC', (SPCSETS(J),J=1,NSUB)
               WRITE(F06,1830) 'SPC', (SPCSETS(J),J=1,NSUB)
            ENDIF
         ENDIF
      ENDDO

      ! Make sure that, if there are more than 1 of MPC or MPC set requests in
      ! Case Control, that all (nonzero) are the same. Otherwise issue error.
      MPCSET = MPCSETS(1)
      DO I=2,NSUB
         IF (MPCSETS(I) /= 0) THEN
            IF (MPCSETS(I) /= MPCSET) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1830) 'MPC', (MPCSETS(J),J=1,NSUB)
               WRITE(F06,1830) 'MPC', (MPCSETS(J),J=1,NSUB)
            ENDIF
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

  901 FORMAT(' *WARNING    : ELDATA ENTRY NOT ALLOWED IN RESTART')

  902 FORMAT(' *WARNING    : CASE CONTROL ENTRY: ',A                                                                               &
                    ,/,14X,' IS NOT ALLOWED. ALL FOLLOWING CASE CONTROL COMMANDS IGNORED')

 1004 FORMAT(' *ERROR  1004: NO METHOD ENTRY FOUND IN CASE CONTROL DECK. CANNOT RUN REAL EIGENVALUE OR LINEAR BUCKLING ANALYSIS.')

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ',A,' ENTRY. ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

 1015 FORMAT(' *WARNING    : GRID POINT FORCE BALANCE ONLY ALLOWED IN ',A,' SOLUTION. ABOVE ENTRY IGNORED')

 1016 FORMAT(' *INFORMATION: ALL SUBCASES WILL USE "',A,'" AS THE LOCATION OF ',A,' OUTPUTS FOR PSHELL QUAD4 ELEMENTS')

 1028 FORMAT(' *ERROR  1028: THERE MUST BE 2 SUBCASES FOR LINEAR BUCKLING ANALYSES BUT NSUB = ',I8)

 1101 FORMAT(' *ERROR  1101: FOR BUCKLING ANALYSES THERE MUST BE AT LEAST ONE SUBCASE WITH A METHOD (BUCKLING EIGENPROBLEM)',     &
                           ' AND AT LEAST ONE OTHER SUBCASE THAT CARRIES A LOAD (AND/OR TEMP) TO ACT AS THE STATSUB PRELOAD.')

 1102 FORMAT(' *ERROR  1102: FOR BUCKLING ANALYSES ANY LOAD, SPS OR MPC IN 2nd SUBCASE MUST BE THE SAME AS THOSE IN 1st SUBCASE')

 1103 FORMAT(' *ERROR  1103: BUCKLING SUBCASE (INTERNAL #',I8,') HAS NO STATSUB AND NO STATIC SUBCASE WITH A LOAD WAS FOUND',     &
                           ' TO ACT AS THE LEGACY DEFAULT PRELOAD.')

 1104 FORMAT(' *ERROR  1104: STATSUB(PRELOAD)=',I8,' ON BUCKLING SUBCASE (INTERNAL #',I8,') REFERENCES A SUBCASE THAT DOES NOT',  &
                           ' EXIST IN THIS CASE CONTROL DECK.')

 1105 FORMAT(' *ERROR  1105: STATSUB(PRELOAD)=',I8,' ON BUCKLING SUBCASE (INTERNAL #',I8,') REFERENCES ANOTHER BUCKLING SUBCASE.',&
                           ' STATSUB MUST POINT TO A LINEAR-STATIC SUBCASE.')

 1106 FORMAT(' *ERROR  1106: STATSUB(PRELOAD)=',I8,' ON BUCKLING SUBCASE (INTERNAL #',I8,') REFERENCES A SUBCASE THAT CARRIES',   &
                           ' NEITHER LOAD NOR TEMP.')

 1199 FORMAT(' *WARNING    : BE CAREFUL WITH LINES THAT BEGIN WITH A $ SIGN IN COL 1 FOLLOWED BY AN UPPER CASE LETTER IN EXEC OR', &
                           ' CASE CONTROL.'                                                                                        &
                    ,/,14X,' THE LINE CAN BE MISINTERPRETED AS A DIRECTIVE FOR THE BANDIT GRID RESEQUENCING ALGORITHM.'            &
                    ,/,14X,' SEE THE BANDIT.PDF FILE INSTALLED WHEN YOU RAN SETUP.EXE TO INSTALL MYSTRAN' &
                    ,/,14X, A)

 1803 FORMAT(' *ERROR  1803: CASE CONTROL FOUND NTSUB = ',I8,' TEMP ENTRIES WHEN THERE ARE ONLY NSUB = ',I8,' SUBCASES.',          &
                           ' NTSUB MUST BE <= NSUB')

 1830 FORMAT(' *ERROR  1830: ONLY ONE ',A,' SET ID IS ALLOWED PER RUN. HOWEVER, IN CASE CONTROL THERE WERE THE FOLLOWING SET',     &
                           ' ID''s FOUND: '                                                                                        &
                           ,/,14X,32767(I8,', '))

 9993 FORMAT(' *WARNING    : PRIOR ENTRY NOT PROCESSED BY ',A)

! **********************************************************************************************************************************

      END SUBROUTINE LOADC
