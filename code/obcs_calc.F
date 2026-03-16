C modified by Jorge Magalhaes jmagalhaes@fc.ul.pt, 2026 Feb. See log file and notes below
C
#include "OBCS_OPTIONS.h"

      SUBROUTINE OBCS_CALC( futureTime, futureIter,
     &                      uVel, vVel, wVel, theta, salt,
     &                      myThid )
C     *==========================================================*
C     | SUBROUTINE OBCS_CALC
C     | o Calculate future boundary data at open boundaries
C     |   at time = futureTime
C     *==========================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_EXCH2
# include "W2_EXCH2_SIZE.h"
#endif /* ALLOW_EXCH2 */
#include "SET_GRID.h"
#include "GRID.h"
#include "OBCS_PARAMS.h"
#include "OBCS_GRID.h"
#include "OBCS_FIELDS.h"
#include "EOS.h"

C     == Routine arguments ==
      INTEGER futureIter
      _RL futureTime
      _RL uVel (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL vVel (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL wVel (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL theta(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL salt (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      INTEGER myThid

#ifdef ALLOW_OBCS

C     == Local variables ==			!JM added 'JJ' here for MB tidal forcing. See *.txt logs
      INTEGER bi, bj
      INTEGER I, J ,K, JJ					
      
C     -- Shared and Tutorial Variables -- 	!JM Added for MB tidal forcing. See *.txt logs
      _RL obTimeScale, Uinflow, rampTime2
      _RL vertStructWst(Nr)
      _RL mz, strat, kx, tmpsum

C     -- Specific Variables for MB 		!JM Added for MB tidal forcing. See *.txt logs
      _RL facdmp, dumvar
      _RL omegaF, ftide
      _RL omgF(8), amplF(8), phsF(8)

C     --------------------------------------------------------------------------------------------------------
C     == Step 2: Tidal Constants and the Anti-Crash Ramp ==	!JM Added for MB tidal forcing. See *.txt logs

      WRITE(standardMessageUnit,*) 'JM runtime check obcs_calc.F CP01'

C     Frequencies (rad/s)
      omgF(1) = 0.000140518900
      omgF(2) = 0.000145444100
      omgF(3) = 0.000137879700
      omgF(4) = 0.000145842300
      omgF(5) = 0.000072921170
      omgF(6) = 0.000067597740
      omgF(7) = 0.000072522950
      omgF(8) = 0.000064958540

C     Amplitudes (m/s)
      amplF(1) = 0.029840145660
      amplF(2) = 0.017675093411
      amplF(3) = 0.004950825738
      amplF(4) = 0.005995579943
      amplF(5) = 0.009220553631
      amplF(6) = 0.005606491338
      amplF(7) = 0.002531745179
      amplF(8) = 0.001195082500

C     Phases (rad), assumes model start time: t=2009     2     2    11    30     0
      phsF(1) = -0.354637125401
      phsF(2) = -2.882378603450
      phsF(3) = -1.498207132120
      phsF(4) = -1.504070929408
      phsF(5) = -2.197846922123
      phsF(6) = 2.575411209660
      phsF(7) = -0.637532956749
      phsF(8) = 1.580049693138

C     The Ramp: Grows from 0 to 1 over rampTime2 seconds
      rampTime2 = 1. _d 4
      facdmp = 1. _d 0 - exp(-1. _d 0/rampTime2*futureTime)

C     --------------------------------------------------------------------------------------------------------

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('OBCS_CALC',myThid)
#endif

      DO bj=myByLo(myThid),myByHi(myThid)
      DO bi=myBxLo(myThid),myBxHi(myThid)

C     Eastern OB
      IF (useOrlanskiEast) THEN
        CALL ORLANSKI_EAST(
     &          bi, bj, futureTime,
     &          uVel, vVel, wVel, theta, salt,
     &          myThid )
      ELSE
        DO K=1,Nr
          DO J=1-OLy,sNy+OLy
	    dumvar = 0. _d 0
            DO JJ=1,8
              dumvar = dumvar
     &        +amplF(JJ)*cos(omgF(JJ)*futureTime-phsF(JJ))*facdmp
            ENDDO
            OBEu(J,K,bi,bj)=dumvar - 0.0141				!JM add mean U here
            OBEv(J,K,bi,bj)=0.
            OBEt(J,K,bi,bj)=tRef(K)
            OBEs(J,K,bi,bj)=sRef(K)
#ifdef ALLOW_NONHYDROSTATIC
            OBEw(J,K,bi,bj)=0.
#endif
          ENDDO
        ENDDO
      ENDIF

C     Western OB
      IF (useOrlanskiWest) THEN
        CALL ORLANSKI_WEST(
     &          bi, bj, futureTime,
     &          uVel, vVel, wVel, theta, salt,
     &          myThid )
      ELSE
        DO K=1,Nr
          DO J=1-OLy,sNy+OLy
            dumvar = 0. _d 0
            DO JJ=1,8
              dumvar = dumvar
     &        +amplF(JJ)*cos(omgF(JJ)*futureTime-phsF(JJ))*facdmp
            ENDDO
            OBWu(J,K,bi,bj)=dumvar - 0.0141				!JM add mean U here
            OBWv(J,K,bi,bj)=0.
            OBWt(J,K,bi,bj)=tRef(K)
            OBWs(J,K,bi,bj)=sRef(K)
#ifdef ALLOW_NONHYDROSTATIC
            OBWw(J,K,bi,bj)=0.
#endif
          ENDDO
        ENDDO
      ENDIF

C         Northern OB, template for forcing
      IF (useOrlanskiNorth) THEN
        CALL ORLANSKI_NORTH(
     &          bi, bj, futureTime,
     &          uVel, vVel, wVel, theta, salt,
     &          myThid )
      ELSE
        DO K=1,Nr
          DO I=1-OLx,sNx+OLx
            OBNv(I,K,bi,bj)=0.
            OBNu(I,K,bi,bj)=0.
            OBNt(I,K,bi,bj)=tRef(K)
            OBNs(I,K,bi,bj)=sRef(K)
#ifdef ALLOW_NONHYDROSTATIC
            OBNw(I,K,bi,bj)=0.
#endif
          ENDDO
        ENDDO
      ENDIF

C         Southern OB, template for forcing
      IF (useOrlanskiSouth) THEN
        CALL ORLANSKI_SOUTH(
     &          bi, bj, futureTime,
     &          uVel, vVel, wVel, theta, salt,
     &          myThid )
      ELSE
        DO K=1,Nr
          DO I=1-OLx,sNx+OLx
            OBSu(I,K,bi,bj)=0.
            OBSv(I,K,bi,bj)=0.
            OBSt(I,K,bi,bj)=tRef(K)
            OBSs(I,K,bi,bj)=sRef(K)
#ifdef ALLOW_NONHYDROSTATIC
            OBSw(I,K,bi,bj)=0.
#endif
          ENDDO
        ENDDO
      ENDIF

C--   end bi,bj loops.
      ENDDO
      ENDDO

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('OBCS_CALC',myThid)
#endif
#endif /* ALLOW_OBCS */

      RETURN
      END
