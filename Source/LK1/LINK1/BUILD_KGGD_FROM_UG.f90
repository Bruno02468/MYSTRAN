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

      SUBROUTINE BUILD_KGGD_FROM_UG

! Re-entrant assembly of the G-set differential stiffness matrix KGGD from the displacement field currently held in UG_COL.
!
! This is the SOL 105 step-2 KGGD assembly path that previously lived inline in LINK1. It has been factored out so that the
! multi-buckling-subcase driver added in Phase 4 can rebuild KGGD once per buckling subcase (each with its own preload UG_COL)
! without re-entering all of LINK1.
!
! Caller responsibilities:
!   * UG_COL must already be loaded with the preload static displacement field for the buckling subcase being assembled
!     (typically by READ_L5A_UG_FOR_SUBCASE). Element ELMDIS calls (gated on OPT(6)=='Y' .AND. LOAD_ISTEP>1) read from UG_COL
!     to form per-element KED contributions.
!   * MPC_IND_GRIDS must remain allocated across repeated calls (SPARSE_KGGD consumes it). The original LINK1 step-2 block
!     deallocated MPC_IND_GRIDS immediately after SPARSE_KGGD; for a single-shot invocation (the legacy path) that dealloc
!     happens in LINK1 just after this routine returns, preserving prior behavior.
!
! Re-entry safety: any pre-existing sparse KGGD (I_KGGD, J_KGGD, KGGD) and STF linked-list arrays (STFKEY, STF3) are deallocated
! before fresh allocation so this routine is safe to call multiple times in a row.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, LTERM_KGGD
      USE PARAMS, ONLY                :  ESP0_PAUSE
      USE SPARSE_MATRICES, ONLY       :  I_KGGD, J_KGGD, KGGD
      USE STF_ARRAYS, ONLY            :  STFKEY, STF3
      USE MODEL_STUF, ONLY            :  AGRID, BGRID

      USE BUILD_KGGD_FROM_UG_USE_IFs
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  LINK_MESSAGE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BUILD_KGGD_FROM_UG'

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)
      CHARACTER( 1*BYTE)              :: RESPONSE          ! Used only if ESP0_PAUSE == 'Y'

      INTEGER(LONG)                   :: LTERM             ! Local copy of LTERM_KGGD for optional interactive override
      LOGICAL                         :: WE_ALLOCD_SINGLE_ELEM_ARRS



! **********************************************************************************************************************************

! 1) Drop any stale sparse KGGD left over from a prior buckling-subcase iteration.

      IF (ALLOCATED(KGGD)   .OR. ALLOCATED(I_KGGD) .OR. ALLOCATED(J_KGGD)) THEN
         CALL DEALLOCATE_SPARSE_MAT ( 'KGGD' )
      ENDIF

! 2) Drop any stale STF linked-list arrays. ALLOCATE_STF_ARRAYS FATALs if its target is already allocated, so deallocate first.

      IF (ALLOCATED(STFKEY)) CALL DEALLOCATE_STF_ARRAYS ( 'STFKEY' )
      IF (ALLOCATED(STF3))   CALL DEALLOCATE_STF_ARRAYS ( 'STF3' )

! 3) Estimate LTERM_KGGD (subr ESP0 sizes the linked-list storage for the element merge pass).

      CALL ESP0
      CALL LINK_MESSAGE('CALCULATE ESTIMATE OF KGGD MATRIX SIZE        ')
      LTERM = LTERM_KGGD

      IF (ESP0_PAUSE == 'Y') THEN
         WRITE(SC1,'(A,A,I12)') ' From ESP0: ', 'LTERM_KGGD', ' = ', LTERM
         WRITE(SC1,'(A,A)') ' Do you want to change ', 'LTERM_KGGD estimate? (Y/N)'
         READ(*,*) RESPONSE
         IF ((RESPONSE == 'Y') .OR. (RESPONSE == 'y')) THEN
            WRITE(SC1,'(A)') 'Enter new LTERM_KGGD'
            WRITE(SC1,*)
            READ(*,*) LTERM
            LTERM_KGGD = LTERM
            WRITE(SC1,'(A,I12)') 'New LTERM_KGGD will be = ', LTERM
         ENDIF
      ENDIF

! 4) Allocate STF linked-list workspace, run ESP (element-by-element KED merge), then condense to sparse KGGD.
!    SINGLE ELEMENT ARRAYS (AGRID, BGRID, DT, etc.) may have already been deallocated by LINK1 if we are being
!    invoked re-entrantly from LINK4 for a multi-buckling-subcase rebuild. Re-allocate them just for this assembly
!    pass and deallocate again on exit, so the legacy single-shot path's allocation state is preserved.

      WE_ALLOCD_SINGLE_ELEM_ARRS = .FALSE.
      IF (.NOT. ALLOCATED(AGRID)) THEN
         CALL ALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS', SUBR_NAME )
         WE_ALLOCD_SINGLE_ELEM_ARRS = .TRUE.
      ENDIF

      CALL LINK_MESSAGE('ALLOCATE MEM FOR STFKEY, STFCOL, STFPNT, STF')
      CALL ALLOCATE_STF_ARRAYS ( 'STFKEY', SUBR_NAME )
      CALL ALLOCATE_STF_ARRAYS ( 'STF3',   SUBR_NAME )

      CALL LINK_MESSAGE('G-SET STIFFNESS MATRIX PROCESSOR            ')
      CALL ESP

      CALL LINK_MESSAGE('SPARSE KGGD PROCESSOR                       ')
      CALL SPARSE_KGGD

      CALL DEALLOCATE_STF_ARRAYS ( 'STFKEY' )
      CALL DEALLOCATE_STF_ARRAYS ( 'STF3' )

      IF (WE_ALLOCD_SINGLE_ELEM_ARRS) THEN
         CALL DEALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS' )
      ENDIF

      WRITE(SC1,*) CR13



      RETURN

      END SUBROUTINE BUILD_KGGD_FROM_UG
