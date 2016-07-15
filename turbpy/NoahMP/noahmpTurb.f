! Turbulence code from Noah MP
SUBROUTINE RAGRB(ITER   ,VAI    ,RHOAIR ,HG     ,TAH    , & !in
                   ZPD    ,Z0MG   ,Z0HG   ,HCAN   ,UC     , & !in
                   Z0H    ,FV     ,CWP    ,VEGTYP ,MPE    , & !in
                   TV     ,MOZG   ,FHG    ,ILOC   ,JLOC   , & !inout
                   RAMG   ,RAHG   ,RAWG   ,RB     )           !out
! --------------------------------------------------------------------------------------------------
! compute under-canopy aerodynamic resistance RAG and leaf boundary layer
! resistance RB
! --------------------------------------------------------------------------------------------------
  USE NOAHMP_VEG_PARAMETERS
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs

  INTEGER,              INTENT(IN) :: ILOC   !grid index
  INTEGER,              INTENT(IN) :: JLOC   !grid index
  INTEGER,              INTENT(IN) :: ITER   !iteration index
  INTEGER,              INTENT(IN) :: VEGTYP !vegetation physiology type
  REAL,                 INTENT(IN) :: VAI    !total LAI + stem area index, one sided
  REAL,                 INTENT(IN) :: RHOAIR !density air (kg/m3)
  REAL,                 INTENT(IN) :: HG     !ground sensible heat flux (w/m2)
  REAL,                 INTENT(IN) :: TV     !vegetation temperature (k)
  REAL,                 INTENT(IN) :: TAH    !air temperature at height z0h+zpd (k)
  REAL,                 INTENT(IN) :: ZPD    !zero plane displacement (m)
  REAL,                 INTENT(IN) :: Z0MG   !roughness length, momentum, ground (m)
  REAL,                 INTENT(IN) :: HCAN   !canopy height (m) [note: hcan >= z0mg]
  REAL,                 INTENT(IN) :: UC     !wind speed at top of canopy (m/s)
  REAL,                 INTENT(IN) :: Z0H    !roughness length, sensible heat (m)
  REAL,                 INTENT(IN) :: Z0HG   !roughness length, sensible heat, ground (m)
  REAL,                 INTENT(IN) :: FV     !friction velocity (m/s)
  REAL,                 INTENT(IN) :: CWP    !canopy wind parameter
  REAL,                 INTENT(IN) :: MPE    !prevents overflow error if division by zero

! in & out

  REAL,              INTENT(INOUT) :: MOZG   !Monin-Obukhov stability parameter
  REAL,              INTENT(INOUT) :: FHG    !stability correction

! outputs
  REAL                             :: RAMG   !aerodynamic resistance for momentum (s/m)
  REAL                             :: RAHG   !aerodynamic resistance for sensible heat (s/m)
  REAL                             :: RAWG   !aerodynamic resistance for water vapor (s/m)
  REAL                             :: RB     !bulk leaf boundary layer resistance (s/m)


  REAL :: KH           !turbulent transfer coefficient, sensible heat, (m2/s)
  REAL :: TMP1         !temporary calculation
  REAL :: TMP2         !temporary calculation
  REAL :: TMPRAH2      !temporary calculation for aerodynamic resistances
  REAL :: TMPRB        !temporary calculation for rb
  real :: MOLG,FHGNEW,CWPC
! --------------------------------------------------------------------------------------------------
! stability correction to below canopy resistance

       MOZG = 0.
       MOLG = 0.

       IF(ITER > 1) THEN
        TMP1 = VKC * (GRAV/TAH) * HG/(RHOAIR*CPAIR)
        IF (ABS(TMP1) .LE. MPE) TMP1 = MPE
        MOLG = -1. * FV**3 / TMP1
        MOZG = MIN( (ZPD-Z0MG)/MOLG, 1.)
       END IF

       IF (MOZG < 0.) THEN
          FHGNEW  = (1. - 15.*MOZG)**(-0.25)
       ELSE
          FHGNEW  = 1.+ 4.7*MOZG
       ENDIF

       IF (ITER == 1) THEN
          FHG = FHGNEW
       ELSE
          FHG = 0.5 * (FHG+FHGNEW)
       ENDIF

       CWPC = (CWP * VAI * HCAN * FHG)**0.5
!       CWPC = (CWP*FHG)**0.5

       TMP1 = EXP( -CWPC*Z0HG/HCAN )
       TMP2 = EXP( -CWPC*(Z0H+ZPD)/HCAN )
       TMPRAH2 = HCAN*EXP(CWPC) / CWPC * (TMP1-TMP2)

! aerodynamic resistances raw and rah between heights zpd+z0h and z0hg.

       KH  = MAX ( VKC*FV*(HCAN-ZPD), MPE )
       RAMG = 0.
       RAHG = TMPRAH2 / KH
       RAWG = RAHG

! leaf boundary layer resistance

       TMPRB  = CWPC*50. / (1. - EXP(-CWPC/2.))
       RB     = TMPRB * SQRT(DLEAF(VEGTYP)/UC)
!       RB = 200

  END SUBROUTINE RAGRB

! ==================================================================================================

  SUBROUTINE SFCDIF1(ITER   ,SFCTMP ,RHOAIR ,H      ,QAIR   , & !in
       &             ZLVL   ,ZPD    ,Z0M    ,Z0H    ,UR     , & !in
       &             MPE    ,ILOC   ,JLOC   ,                 & !in
       &             MOZ    ,MOZSGN ,FM     ,FH     ,FM2,FH2, & !inout
       &             CM     ,CH     ,FV     ,CH2     )          !out
! -------------------------------------------------------------------------------------------------
! computing surface drag coefficient CM for momentum and CH for heat
! -------------------------------------------------------------------------------------------------
    IMPLICIT NONE
! -------------------------------------------------------------------------------------------------
! inputs
    
    INTEGER,              INTENT(IN) :: ILOC   !grid index
    INTEGER,              INTENT(IN) :: JLOC   !grid index
    INTEGER,              INTENT(IN) :: ITER   !iteration index
    REAL,                 INTENT(IN) :: SFCTMP !temperature at reference height (k)
    REAL,                 INTENT(IN) :: RHOAIR !density air (kg/m**3)
    REAL,                 INTENT(IN) :: H      !sensible heat flux (w/m2) [+ to atm]
    REAL,                 INTENT(IN) :: QAIR   !specific humidity at reference height (kg/kg)
    REAL,                 INTENT(IN) :: ZLVL   !reference height  (m)
    REAL,                 INTENT(IN) :: ZPD    !zero plane displacement (m)
    REAL,                 INTENT(IN) :: Z0H    !roughness length, sensible heat, ground (m)
    REAL,                 INTENT(IN) :: Z0M    !roughness length, momentum, ground (m)
    REAL,                 INTENT(IN) :: UR     !wind speed (m/s)
    REAL,                 INTENT(IN) :: MPE    !prevents overflow error if division by zero
! in & out

    INTEGER,           INTENT(INOUT) :: MOZSGN !number of times moz changes sign
    REAL,              INTENT(INOUT) :: MOZ    !Monin-Obukhov stability (z/L)
    REAL,              INTENT(INOUT) :: FM     !momentum stability correction, weighted by prior iters
    REAL,              INTENT(INOUT) :: FH     !sen heat stability correction, weighted by prior iters
    REAL,              INTENT(INOUT) :: FM2    !sen heat stability correction, weighted by prior iters
    REAL,              INTENT(INOUT) :: FH2    !sen heat stability correction, weighted by prior iters

! outputs

    REAL,                INTENT(OUT) :: CM     !drag coefficient for momentum
    REAL,                INTENT(OUT) :: CH     !drag coefficient for heat
    REAL,                INTENT(OUT) :: FV     !friction velocity (m/s)
    REAL,                INTENT(OUT) :: CH2    !drag coefficient for heat

! locals
    REAL    :: MOL                      !Monin-Obukhov length (m)
    REAL    :: TMPCM                    !temporary calculation for CM
    REAL    :: TMPCH                    !temporary calculation for CH
    REAL    :: FMNEW                    !stability correction factor, momentum, for current moz
    REAL    :: FHNEW                    !stability correction factor, sen heat, for current moz
    REAL    :: MOZOLD                   !Monin-Obukhov stability parameter from prior iteration
    REAL    :: TMP1,TMP2,TMP3,TMP4,TMP5 !temporary calculation
    REAL    :: TVIR                     !temporary virtual temperature (k)
    REAL    :: MOZ2                     !2/L
    REAL    :: TMPCM2                   !temporary calculation for CM2
    REAL    :: TMPCH2                   !temporary calculation for CH2
    REAL    :: FM2NEW                   !stability correction factor, momentum, for current moz
    REAL    :: FH2NEW                   !stability correction factor, sen heat, for current moz
    REAL    :: TMP12,TMP22,TMP32        !temporary calculation

    REAL    :: CMFM, CHFH, CM2FM2, CH2FH2
! -------------------------------------------------------------------------------------------------
! Monin-Obukhov stability parameter moz for next iteration

    MOZOLD = MOZ
  
    IF(ZLVL <= ZPD) THEN
       write(*,*) 'critical problem: ZLVL <= ZPD; model stops'
       call wrf_error_fatal("STOP in Noah-MP")
    ENDIF

    TMPCM = LOG((ZLVL-ZPD) / Z0M)
    TMPCH = LOG((ZLVL-ZPD) / Z0H)
    TMPCM2 = LOG((2.0 + Z0M) / Z0M)
    TMPCH2 = LOG((2.0 + Z0H) / Z0H)

    IF(ITER == 1) THEN
       FV   = 0.0
       MOZ  = 0.0
       MOL  = 0.0
       MOZ2 = 0.0
    ELSE
       TVIR = (1. + 0.61*QAIR) * SFCTMP
       TMP1 = VKC * (GRAV/TVIR) * H/(RHOAIR*CPAIR)
       IF (ABS(TMP1) .LE. MPE) TMP1 = MPE
       MOL  = -1. * FV**3 / TMP1
       MOZ  = MIN( (ZLVL-ZPD)/MOL, 1.)
       MOZ2  = MIN( (2.0 + Z0H)/MOL, 1.)
    ENDIF

! accumulate number of times moz changes sign.

    IF (MOZOLD*MOZ .LT. 0.) MOZSGN = MOZSGN+1
    IF (MOZSGN .GE. 2) THEN
       MOZ = 0.
       FM = 0.
       FH = 0.
       MOZ2 = 0.
       FM2 = 0.
       FH2 = 0.
    ENDIF

! evaluate stability-dependent variables using moz from prior iteration
    IF (MOZ .LT. 0.) THEN
       TMP1 = (1. - 16.*MOZ)**0.25
       TMP2 = LOG((1.+TMP1*TMP1)/2.)
       TMP3 = LOG((1.+TMP1)/2.)
       FMNEW = 2.*TMP3 + TMP2 - 2.*ATAN(TMP1) + 1.5707963
       FHNEW = 2*TMP2

! 2-meter
       TMP12 = (1. - 16.*MOZ2)**0.25
       TMP22 = LOG((1.+TMP12*TMP12)/2.)
       TMP32 = LOG((1.+TMP12)/2.)
       FM2NEW = 2.*TMP32 + TMP22 - 2.*ATAN(TMP12) + 1.5707963
       FH2NEW = 2*TMP22
    ELSE
       FMNEW = -5.*MOZ
       FHNEW = FMNEW
       FM2NEW = -5.*MOZ2
       FH2NEW = FM2NEW
    ENDIF

! except for first iteration, weight stability factors for previous
! iteration to help avoid flip-flops from one iteration to the next

    IF (ITER == 1) THEN
       FM = FMNEW
       FH = FHNEW
       FM2 = FM2NEW
       FH2 = FH2NEW
    ELSE
       FM = 0.5 * (FM+FMNEW)
       FH = 0.5 * (FH+FHNEW)
       FM2 = 0.5 * (FM2+FM2NEW)
       FH2 = 0.5 * (FH2+FH2NEW)
    ENDIF

! exchange coefficients

    FH = MIN(FH,0.9*TMPCH)
    FM = MIN(FM,0.9*TMPCM)
    FH2 = MIN(FH2,0.9*TMPCH2)
    FM2 = MIN(FM2,0.9*TMPCM2)

    CMFM = TMPCM-FM
    CHFH = TMPCH-FH
    CM2FM2 = TMPCM2-FM2
    CH2FH2 = TMPCH2-FH2
    IF(ABS(CMFM) <= MPE) CMFM = MPE
    IF(ABS(CHFH) <= MPE) CHFH = MPE
    IF(ABS(CM2FM2) <= MPE) CM2FM2 = MPE
    IF(ABS(CH2FH2) <= MPE) CH2FH2 = MPE
    CM  = VKC*VKC/(CMFM*CMFM)
    CH  = VKC*VKC/(CMFM*CHFH)
    CH2  = VKC*VKC/(CM2FM2*CH2FH2)
        
! friction velocity

    FV = UR * SQRT(CM)
    CH2  = VKC*FV/CH2FH2

  END SUBROUTINE SFCDIF1

! ==================================================================================================

  SUBROUTINE SFCDIF2(ITER   ,Z0     ,THZ0   ,THLM   ,SFCSPD , & !in
                     CZIL   ,ZLM    ,ILOC   ,JLOC   ,         & !in
                     AKMS   ,AKHS   ,RLMO   ,WSTAR2 ,         & !in
                     USTAR  )                                   !out

! -------------------------------------------------------------------------------------------------
! SUBROUTINE SFCDIF (renamed SFCDIF_off to avoid clash with Eta PBL)
! -------------------------------------------------------------------------------------------------
! CALCULATE SURFACE LAYER EXCHANGE COEFFICIENTS VIA ITERATIVE PROCESS.
! SEE CHEN ET AL (1997, BLM)
! -------------------------------------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ILOC
    INTEGER, INTENT(IN) :: JLOC
    INTEGER, INTENT(IN) :: ITER
    REAL,    INTENT(IN) :: ZLM, Z0, THZ0, THLM, SFCSPD, CZIL
    REAL, intent(INOUT) :: AKMS
    REAL, intent(INOUT) :: AKHS
    REAL, intent(INOUT) :: RLMO
    REAL, intent(INOUT) :: WSTAR2
    REAL,   intent(OUT) :: USTAR

    REAL     ZZ, PSLMU, PSLMS, PSLHU, PSLHS
    REAL     XX, PSPMU, YY, PSPMS, PSPHU, PSPHS
    REAL     ZILFC, ZU, ZT, RDZ, CXCH
    REAL     DTHV, DU2, BTGH, ZSLU, ZSLT, RLOGU, RLOGT
    REAL     ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4

    REAL     XLU, XLT, XU, XT, PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN,  &
         &         RLMA

    INTEGER  ILECH, ITR

    INTEGER, PARAMETER :: ITRMX  = 5
    REAL,    PARAMETER :: WWST   = 1.2
    REAL,    PARAMETER :: WWST2  = WWST * WWST
    REAL,    PARAMETER :: VKRM   = 0.40
    REAL,    PARAMETER :: EXCM   = 0.001
    REAL,    PARAMETER :: BETA   = 1.0 / 270.0
    REAL,    PARAMETER :: BTG    = BETA * GRAV
    REAL,    PARAMETER :: ELFC   = VKRM * BTG
    REAL,    PARAMETER :: WOLD   = 0.15
    REAL,    PARAMETER :: WNEW   = 1.0 - WOLD
    REAL,    PARAMETER :: PIHF   = 3.14159265 / 2.
    REAL,    PARAMETER :: EPSU2  = 1.E-4
    REAL,    PARAMETER :: EPSUST = 0.07
    REAL,    PARAMETER :: EPSIT  = 1.E-4
    REAL,    PARAMETER :: EPSA   = 1.E-8
    REAL,    PARAMETER :: ZTMIN  = -5.0
    REAL,    PARAMETER :: ZTMAX  = 1.0
    REAL,    PARAMETER :: HPBL   = 1000.0
    REAL,    PARAMETER :: SQVISC = 258.2
    REAL,    PARAMETER :: RIC    = 0.183
    REAL,    PARAMETER :: RRIC   = 1.0 / RIC
    REAL,    PARAMETER :: FHNEU  = 0.8
    REAL,    PARAMETER :: RFC    = 0.191
    REAL,    PARAMETER :: RFAC   = RIC / ( FHNEU * RFC * RFC )

! ----------------------------------------------------------------------
! NOTE: THE TWO CODE BLOCKS BELOW DEFINE FUNCTIONS
! ----------------------------------------------------------------------
! LECH'S SURFACE FUNCTIONS
    PSLMU (ZZ)= -0.96* log (1.0-4.5* ZZ)
    PSLMS (ZZ)= ZZ * RRIC -2.076* (1. -1./ (ZZ +1.))
    PSLHU (ZZ)= -0.96* log (1.0-4.5* ZZ)
    PSLHS (ZZ)= ZZ * RFAC -2.076* (1. -1./ (ZZ +1.))
! PAULSON'S SURFACE FUNCTIONS
    PSPMU (XX)= -2.* log ( (XX +1.)*0.5) - log ( (XX * XX +1.)*0.5)   &
         &        +2.* ATAN (XX)                                            &
         &- PIHF
    PSPMS (YY)= 5.* YY
    PSPHU (XX)= -2.* log ( (XX * XX +1.)*0.5)
    PSPHS (YY)= 5.* YY

! THIS ROUTINE SFCDIF CAN HANDLE BOTH OVER OPEN WATER (SEA, OCEAN) AND
! OVER SOLID SURFACE (LAND, SEA-ICE).
! ----------------------------------------------------------------------
!     ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
!     C......ZTFC=0.1
!     CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
! ----------------------------------------------------------------------
    ILECH = 0

! ----------------------------------------------------------------------
    ZILFC = - CZIL * VKRM * SQVISC
    ZU = Z0
    RDZ = 1./ ZLM
    CXCH = EXCM * RDZ
    DTHV = THLM - THZ0

! BELJARS CORRECTION OF USTAR
    DU2 = MAX (SFCSPD * SFCSPD,EPSU2)
    BTGH = BTG * HPBL

    IF(ITER == 1) THEN
        IF (BTGH * AKHS * DTHV .ne. 0.0) THEN
           WSTAR2 = WWST2* ABS (BTGH * AKHS * DTHV)** (2./3.)
        ELSE
           WSTAR2 = 0.0
        END IF
        USTAR = MAX (SQRT (AKMS * SQRT (DU2+ WSTAR2)),EPSUST)
        RLMO = ELFC * AKHS * DTHV / USTAR **3
    END IF
 
! ZILITINKEVITCH APPROACH FOR ZT
    ZT = MAX(1.E-6,EXP (ZILFC * SQRT (USTAR * Z0))* Z0)
    ZSLU = ZLM + ZU
    ZSLT = ZLM + ZT
    RLOGU = log (ZSLU / ZU)
    RLOGT = log (ZSLT / ZT)

! ----------------------------------------------------------------------
! 1./MONIN-OBUKKHOV LENGTH-SCALE
! ----------------------------------------------------------------------
    ZETALT = MAX (ZSLT * RLMO,ZTMIN)
    RLMO = ZETALT / ZSLT
    ZETALU = ZSLU * RLMO
    ZETAU = ZU * RLMO
    ZETAT = ZT * RLMO

    IF (ILECH .eq. 0) THEN
       IF (RLMO .lt. 0.)THEN
          XLU4 = 1. -16.* ZETALU
          XLT4 = 1. -16.* ZETALT
          XU4  = 1. -16.* ZETAU
          XT4  = 1. -16.* ZETAT
          XLU  = SQRT (SQRT (XLU4))
          XLT  = SQRT (SQRT (XLT4))
          XU   = SQRT (SQRT (XU4))

          XT = SQRT (SQRT (XT4))
          PSMZ = PSPMU (XU)
          SIMM = PSPMU (XLU) - PSMZ + RLOGU
          PSHZ = PSPHU (XT)
          SIMH = PSPHU (XLT) - PSHZ + RLOGT
       ELSE
          ZETALU = MIN (ZETALU,ZTMAX)
          ZETALT = MIN (ZETALT,ZTMAX)
          PSMZ = PSPMS (ZETAU)
          SIMM = PSPMS (ZETALU) - PSMZ + RLOGU
          PSHZ = PSPHS (ZETAT)
          SIMH = PSPHS (ZETALT) - PSHZ + RLOGT
       END IF
! ----------------------------------------------------------------------
! LECH'S FUNCTIONS
! ----------------------------------------------------------------------
    ELSE
       IF (RLMO .lt. 0.)THEN
          PSMZ = PSLMU (ZETAU)
          SIMM = PSLMU (ZETALU) - PSMZ + RLOGU
          PSHZ = PSLHU (ZETAT)
          SIMH = PSLHU (ZETALT) - PSHZ + RLOGT
       ELSE
          ZETALU = MIN (ZETALU,ZTMAX)
          ZETALT = MIN (ZETALT,ZTMAX)
          PSMZ = PSLMS (ZETAU)
          SIMM = PSLMS (ZETALU) - PSMZ + RLOGU
          PSHZ = PSLHS (ZETAT)
          SIMH = PSLHS (ZETALT) - PSHZ + RLOGT
       END IF
! ----------------------------------------------------------------------
       END IF

! ----------------------------------------------------------------------
! BELJAARS CORRECTION FOR USTAR
! ----------------------------------------------------------------------
       USTAR = MAX (SQRT (AKMS * SQRT (DU2+ WSTAR2)),EPSUST)

! ZILITINKEVITCH FIX FOR ZT
       ZT = MAX(1.E-6,EXP (ZILFC * SQRT (USTAR * Z0))* Z0)
       ZSLT = ZLM + ZT
!-----------------------------------------------------------------------
       RLOGT = log (ZSLT / ZT)
       USTARK = USTAR * VKRM
       AKMS = MAX (USTARK / SIMM,CXCH)
!-----------------------------------------------------------------------
! IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
!-----------------------------------------------------------------------
       AKHS = MAX (USTARK / SIMH,CXCH)

       IF (BTGH * AKHS * DTHV .ne. 0.0) THEN
          WSTAR2 = WWST2* ABS (BTGH * AKHS * DTHV)** (2./3.)
       ELSE
          WSTAR2 = 0.0
       END IF
!-----------------------------------------------------------------------
       RLMN = ELFC * AKHS * DTHV / USTAR **3
!-----------------------------------------------------------------------
!     IF(ABS((RLMN-RLMO)/RLMA).LT.EPSIT)    GO TO 110
!-----------------------------------------------------------------------
       RLMA = RLMO * WOLD+ RLMN * WNEW
!-----------------------------------------------------------------------
       RLMO = RLMA

!       write(*,'(a20,10f15.6)')'SFCDIF: RLMO=',RLMO,RLMN,ELFC , AKHS , DTHV , USTAR
!    END DO
! ----------------------------------------------------------------------
  END SUBROUTINE SFCDIF2

! ==================================================================================================
  SUBROUTINE SFCDIF3(ILOC   ,JLOC   ,TSK    ,QS     ,PSFC   ,&  !in  
                     PBLH   ,Z0     ,Z0BASE ,VEGTYP ,ISURBAN,&  !in  
                     IZ0TLND,SFCSPD ,ITER   ,ITRMX  ,TLOW   ,&  !in  
                     THLOW  ,QLOW   ,CWMLOW ,ZSL    ,        &  !in
                     PLOW   ,USTAR  ,AKMS   ,AKHS   ,CHS2   ,&  !inout
                     CQS2   ,RLMO   )                           !out 

  USE MODULE_SF_MYJSFC, ONLY :   &
     &                         EPSU2  , &  
     &                         EPSUST , &
     &                         EPSZT  , &  
     &                         BETA   , &  
     &                         EXCML  , &  
     &                         RIC    , &  
     &                         SQVISC , &
     &                         ZTFC   , &  
     &                         BTG    , &  
     &                         CZIV   , &  
     &                         PI     , &  
     &                         PIHF   , &  
     &                         KZTM   , &  
     &                         KZTM2  , &  
     &                         DZETA1 , &
     &                         DZETA2 , &
     &                         FH01   , &  
     &                         FH02   , &  
     &                         WWST2  , &  
     &                         WWST   , &  
     &                         ZTMAX1 , &
     &                         ZTMAX2 , &
     &                         ZTMIN1 , &
     &                         ZTMIN2 , &
     &                         PSIH1  , &  
     &                         PSIH2  , &  
     &                         PSIM1  , &  
     &                         PSIM2   

  USE MODULE_MODEL_CONSTANTS

!----------------------------------------------------------------------   
!  computing surface drag coefficient CM for momentum and CH for heat
!  Joakim Refslund, 2011, MYJ SFCLAY
!----------------------------------------------------------------------   
   IMPLICIT NONE     
!----------------------------------------------------------------------   
! input
    INTEGER,INTENT(IN) :: ILOC
    INTEGER,INTENT(IN) :: JLOC
    REAL   ,INTENT(IN) :: TSK
    REAL   ,INTENT(IN) :: PSFC
    REAL   ,INTENT(IN) :: PBLH
    INTEGER,INTENT(IN) :: VEGTYP  !in routine
    INTEGER,INTENT(IN) :: ISURBAN !in veg_parm                    
    INTEGER,INTENT(IN) :: IZ0TLND
    REAL   ,INTENT(IN) :: QLOW
    REAL   ,INTENT(IN) :: THLOW
    REAL   ,INTENT(IN) :: TLOW
    REAL   ,INTENT(IN) :: CWMLOW
    REAL   ,INTENT(IN) :: SFCSPD
    REAL   ,INTENT(IN) :: PLOW
    REAL   ,INTENT(IN) :: ZSL
    REAL   ,INTENT(IN) :: Z0BASE
    INTEGER,INTENT(IN) :: ITER
    INTEGER,INTENT(IN) :: ITRMX

! output                                                                        
    REAL   ,INTENT(OUT) :: CHS2
    REAL   ,INTENT(OUT) :: CQS2
    REAL   ,INTENT(OUT) :: RLMO

! input/output
    REAL   ,INTENT(INOUT) :: AKHS
    REAL   ,INTENT(INOUT) :: AKMS
    REAL   :: QZ0
    REAL   ,INTENT(INOUT) :: USTAR
    REAL   ,INTENT(IN) :: Z0
    REAL   ,INTENT(INOUT):: QS
    REAL                :: RIB

! local
    INTEGER :: ITR,K
    REAL :: THZ0
    REAL :: THVLOW
    REAL :: CT
    REAL :: BTGH
    REAL :: BTGX
    REAL :: CXCHL
    REAL :: DTHV
    REAL :: DU2
    REAL :: ELFC
    REAL :: PSH02
    REAL :: PSH10
    REAL :: PSHZ
    REAL :: PSHZL
    REAL :: PSM10
    REAL :: PSMZ
    REAL :: PSMZL
    REAL :: RDZ
    REAL :: RDZT
    REAL :: RLMA !???
    REAL :: RLMN !???
    REAL :: RLOGT
    REAL :: RLOGU
    REAL :: RZ
    REAL :: SIMH
    REAL :: SIMM
    REAL :: USTARK
    REAL :: WSTAR2
    REAL :: WSTAR
    REAL :: CHS
    REAL :: RZSU
    REAL :: RZST
    REAL :: X,XLT,XLT4,XLU,XLU4,XT,XT4,XU,XU4,ZETALT,ZETALU         , &
            ZETAT,ZETAU,ZQ,ZSLT,ZSLU,ZT,ZU,TOPOTERM,ZZIL
    REAL :: AKHS02,AKHS10,AKMS02,AKMS10
    REAL :: ZU10
    REAL :: ZT02
    REAL :: ZT10
    REAL :: RLNU10
    REAL :: RLNT02
    REAL :: RLNT10
    REAL :: ZTAU10
    REAL :: ZTAT02
    REAL :: ZTAT10
    REAL :: SIMM10
    REAL :: SIMH02
    REAL :: SIMH10
    REAL :: ZUUZ
    REAL :: EKMS10
    REAL :: test
    REAL :: E1

   REAL,    PARAMETER :: VKRM    = 0.40
   REAL,    PARAMETER :: CZETMAX = 10.

! diagnostic terms                                                         

   REAL :: CZIL
   REAL :: ZILFC

! KTMZ,KTMZ2,DZETA1,DZETA2,FH01,FH02,ZTMAX1,ZTMAX2,ZTMIN1,ZTMIN2,
! PSIH1,PSIH2,PSIM1,PSIM2 ARE DEFINED IN MODULE_SF_MYJSFC

!----------------------------------------------------------------------   
!    IF (ILOC.eq.39 .and. JLOC.eq.63 .and. ITER == 1 ) then
!      write(*,*) "THZ0=",THZ0
!      write(*,*) "QS =",QS
!      write(*,*) "PSFC=",PSFC
!      write(*,*) "PBLH=",PBLH
!      write(*,*) "Z0=",Z0
!      write(*,*) "Z0BASE=",Z0BASE
!      write(*,*) "VEGTYP=",VEGTYP
!      write(*,*) "ISURBAN=",ISURBAN
!      write(*,*) "IZ0TLND=",IZ0TLND
!      write(*,*) "SFCSPD=",SFCSPD
!      write(*,*) "TLOW=",TLOW
!      write(*,*) "THLOW=",THLOW
!      write(*,*) "THVLOW=",THVLOW
!      write(*,*) "QLOW=",QLOW
!      write(*,*) "CWMLOW=",CWMLOW
!      write(*,*) "ZSL=",ZSL
!      write(*,*) "PLOW=",PLOW
!      write(*,*) "USTAR=",USTAR
!      write(*,*) "AKMS=",AKMS
!      write(*,*) "AKHS=",AKHS
!      write(*,*) "CHS2=",CHS2
!      write(*,*) "CQS2=",CQS2
!      write(*,*) "RLMO=",RLMO
!      write(*,*) "ITER=",ITER
!      call wrf_error_fatal("STOP in SFCDIF3")
!    ENDIF

! calculate potential and virtual potential temperatures
    THVLOW = THLOW*(1.+EP_1*QLOW)
    THZ0   = TSK*(P1000mb/PSFC)**RCP

! calculate initial values
    ZU = Z0
    ZT = ZU*ZTFC      !ZTFC = ZOH/ZOM =<1 set to 1 at beginning                                                 
    ZQ = ZT
    QZ0 = QS

    RDZ = 1./ZSL
    CXCHL = EXCML*RDZ
    DTHV = THVLOW-THZ0*(0.608*QZ0+1.)    !delta pot. virtual temperature 

    BTGX=GRAV/THLOW
    ELFC=VKRM*BTGX

! Minimum PBLH is >= 1000.
    IF(PBLH > 1000.)THEN  
       BTGH = BTGX*PBLH
    ELSE
       BTGH = BTGX*1000.
    ENDIF

    DU2 = MAX(SFCSPD*SFCSPD,EPSU2)  !Wind speed - EPSU2 parm = 1*10^-6                    
    RIB = BTGX*DTHV*ZSL/DU2         !Bulk richardson stability                               

    ZSLU = ZSL+ZU
    RZSU = ZSLU/ZU
    RLOGU = LOG(RZSU)       !log(z/z0) 

    ZSLT = ZSL + ZU

    IF ( (IZ0TLND==0) .or. (VEGTYP == ISURBAN) ) THEN    ! ARE IZ0TLND DEFINED HERE?           
        ! Just use the original CZIL value.                          
        CZIL = 0.1
    ELSE
        ! Modify CZIL according to Chen & Zhang, 2009                
        ! CZIL = 10 ** -0.40 H, ( where H = 10*Zo )                  
        CZIL = 10.0 ** ( -0.40 * ( Z0 / 0.07 ) )
    ENDIF
    ZILFC=-CZIL*VKRM*SQVISC     !SQVISC parm                                 

! stable                                                                  
    IF(DTHV>0.)THEN
       IF (RIB<RIC) THEN
          ZZIL=ZILFC*(1.0+(RIB/RIC)*(RIB/RIC)*CZETMAX)
       ELSE
          ZZIL=ZILFC*(1.0+CZETMAX)
       ENDIF
! unstable                                                                
    ELSE
       ZZIL=ZILFC
    ENDIF

!---  ZILITINKEVITCH FIX FOR ZT                                           
! oldform   ZT=MAX(EXP(ZZIL*SQRT(USTAR*ZU))*ZU,EPSZT)                     
    ZT=MAX(EXP(ZZIL*SQRT(USTAR*Z0BASE))*Z0BASE,EPSZT)  !Z0 is backgrund roughness?           
    RZST=ZSLT/ZT
    RLOGT=LOG(RZST)

!----------------------------------------------------------------------   
!  1./MONIN-OBUKHOV LENGTH-SCALE                                       
!----------------------------------------------------------------------   
    RLMO=ELFC*AKHS*DTHV/USTAR**3

    ZETALU=ZSLU*RLMO
    ZETALT=ZSLT*RLMO
    ZETAU=ZU*RLMO
    ZETAT=ZT*RLMO

    ZETALU=MIN(MAX(ZETALU,ZTMIN2),ZTMAX2)
    ZETALT=MIN(MAX(ZETALT,ZTMIN2),ZTMAX2)
    ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)
    ZETAT=MIN(MAX(ZETAT,ZTMIN2/RZST),ZTMAX2/RZST)

!----------------------------------------------------------------------   
!***  LAND FUNCTIONS                                                      
!----------------------------------------------------------------------   

    RZ=(ZETAU-ZTMIN2)/DZETA2
    K=INT(RZ)
    RDZT=RZ-REAL(K)
    K=MIN(K,KZTM2)
    K=MAX(K,0)
    PSMZ=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)

    RZ=(ZETALU-ZTMIN2)/DZETA2
    K=INT(RZ)
    RDZT=RZ-REAL(K)
    K=MIN(K,KZTM2)
    K=MAX(K,0)
    PSMZL=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)

    SIMM=PSMZL-PSMZ+RLOGU

    RZ=(ZETAT-ZTMIN2)/DZETA2
    K=INT(RZ)
    RDZT=RZ-REAL(K)
    K=MIN(K,KZTM2)
    K=MAX(K,0)
    PSHZ=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)

    RZ=(ZETALT-ZTMIN2)/DZETA2
    K=INT(RZ)
    RDZT=RZ-REAL(K)
    K=MIN(K,KZTM2)
    K=MAX(K,0)
    PSHZL=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)

    SIMH=(PSHZL-PSHZ+RLOGT)*FH02
!----------------------------------------------------------------------   
    USTARK=USTAR*VKRM
    AKMS=MAX(USTARK/SIMM,CXCHL)
    AKHS=MAX(USTARK/SIMH,CXCHL)

!----------------------------------------------------------------------   
!  BELJAARS CORRECTION FOR USTAR                                       
!----------------------------------------------------------------------   

    IF(DTHV<=0.)THEN                                           !zj
      WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)                !zj
    ELSE                                                       !zj
      WSTAR2=0.                                                !zj
    ENDIF                                                      !zj

    USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)

    CT=0.
!----------------------------------------------------------------------   
!***  THE FOLLOWING DIAGNOSTIC BLOCK PRODUCES 2-m and 10-m VALUES         
!***  FOR TEMPERATURE, MOISTURE, AND WINDS.  IT IS DONE HERE SINCE        
!***  THE VARIOUS QUANTITIES NEEDED FOR THE COMPUTATION ARE LOST          
!***  UPON EXIT FROM THE ROTUINE.                                         
!----------------------------------------------------------------------   
    WSTAR=SQRT(WSTAR2)/WWST

!jref: calculate in last iteration
!  IF (ITER == ITRMX) THEN

    ZU10=ZU+10.
    ZT02=ZT+02.
    ZT10=ZT+10.

    RLNU10=LOG(ZU10/ZU)
    RLNT02=LOG(ZT02/ZT)
    RLNT10=LOG(ZT10/ZT)

    ZTAU10=ZU10*RLMO
    ZTAT02=ZT02*RLMO
    ZTAT10=ZT10*RLMO

    ZTAU10=MIN(MAX(ZTAU10,ZTMIN2),ZTMAX2)
    ZTAT02=MIN(MAX(ZTAT02,ZTMIN2),ZTMAX2)
    ZTAT10=MIN(MAX(ZTAT10,ZTMIN2),ZTMAX2)

!jref: land diagnostic functions
    RZ=(ZTAU10-ZTMIN2)/DZETA2
    K=INT(RZ)
    RDZT=RZ-REAL(K)
    K=MIN(K,KZTM2)
    K=MAX(K,0)
    PSM10=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)

    SIMM10=PSM10-PSMZ+RLNU10

    RZ=(ZTAT02-ZTMIN2)/DZETA2
    K=INT(RZ)
    RDZT=RZ-REAL(K)
    K=MIN(K,KZTM2)
    K=MAX(K,0)
    PSH02=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)

    SIMH02=(PSH02-PSHZ+RLNT02)*FH02

    RZ=(ZTAT10-ZTMIN2)/DZETA2
    K=INT(RZ)
    RDZT=RZ-REAL(K)
    K=MIN(K,KZTM2)
    K=MAX(K,0)
    PSH10=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)

    SIMH10=(PSH10-PSHZ+RLNT10)*FH02

!jref: diagnostic exchange coefficients                                                                    
    AKMS10=MAX(USTARK/SIMM10,CXCHL)
    AKHS02=MAX(USTARK/SIMH02,CXCHL)
    AKHS10=MAX(USTARK/SIMH10,CXCHL)

!jref: calculation of diagnostics for wind, humidity

!    WSTAR=SQRT(WSTAR2)/WWST
!
!    UMFLX=AKMS*(ULOW -UZ0 )
!    VMFLX=AKMS*(VLOW -VZ0 )
!    HSFLX=AKHS*(THLOW-THZ0)
!    HLFLX=AKHS*(QLOW -QZ0 )

!uncommented for now...
!    U10 =UMFLX/AKMS10+UZ0                                               
!    V10 =VMFLX/AKMS10+VZ0                                               
!    TH02=HSFLX/AKHS02+THZ0                                              
!                                                                         
!***  BE CERTAIN THAT THE 2-M THETA AND 10-M THETA ARE BRACKETED BY       
!***  THE VALUES OF THZ0 AND THLOW.                                       
!                                                                         
!      IF(THLOW>THZ0.AND.(TH02<THZ0.OR.TH02>THLOW).OR.                  &  
!         THLOW<THZ0.AND.(TH02>THZ0.OR.TH02<THLOW))THEN                    
!           TH02=THZ0+2.*RDZ*(THLOW-THZ0)                                  
!      ENDIF                                                               
!                                                                         
!uncommented for now
!      TH10=HSFLX/AKHS10+THZ0                                              
!                                                                         
!      IF(THLOW>THZ0.AND.(TH10<THZ0.OR.TH10>THLOW).OR.                  &  
!         THLOW<THZ0.AND.(TH10>THZ0.OR.TH10<THLOW))THEN                    
!           TH10=THZ0+10.*RDZ*(THLOW-THZ0)                                 
!      ENDIF                                                               
!                                                                         
!      Q02 =HLFLX/AKHS02+QZ0                                               
!      Q10 =HLFLX/AKHS10+QZ0                                               
!jref commented out     
!      TERM1=-0.068283/TLOW                                                
!      PSHLTR=PSFC*EXP(TERM1)                                              
!                                
!----------------------------------------------------------------------   
!***  COMPUTE "EQUIVALENT" Z0 TO APPROXIMATE LOCAL SHELTER READINGS.      
!----------------------------------------------------------------------   
!                                                                         
!      U10E=U10                                                            
!      V10E=V10                                                            
!                                                                         
!      IF(SEAMASK<0.5)THEN                                                 
!LAND :                                                                          
!1st        ZUUZ=MIN(0.5*ZU,0.1)                                          
!1st        ZU=MAX(0.1*ZU,ZUUZ)                                           
!tst        ZUUZ=amin1(ZU*0.50,0.3)                                       
!tst        ZU=amax1(ZU*0.3,ZUUZ)                                         

        ZUUZ=AMIN1(ZU*0.50,0.18)
        ZU=AMAX1(ZU*0.35,ZUUZ)
!                                                                         
        ZU10=ZU+10.
        RZSU=ZU10/ZU
        RLNU10=LOG(RZSU)

        ZETAU=ZU*RLMO
        ZTAU10=ZU10*RLMO

        ZTAU10=MIN(MAX(ZTAU10,ZTMIN2),ZTMAX2)
        ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)

        RZ=(ZTAU10-ZTMIN2)/DZETA2
        K=INT(RZ)
        RDZT=RZ-REAL(K)
        K=MIN(K,KZTM2)
        K=MAX(K,0)
        PSM10=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
        SIMM10=PSM10-PSMZ+RLNU10
        EKMS10=MAX(USTARK/SIMM10,CXCHL)

!        U10E=UMFLX/EKMS10+UZ0                                             
!        V10E=VMFLX/EKMS10+VZ0                                             

!      ENDIF                                                               
!                                                                         
!      U10=U10E                                                            
!      V10=V10E                                                            
!                                                                         
!----------------------------------------------------------------------   
!***  SET OTHER WRF DRIVER ARRAYS                                         
!----------------------------------------------------------------------   
!                                                                         
!jref commented out
!      RLOW=PLOW/(R_D*TLOW)                                                
      CHS=AKHS
      CHS2=AKHS02
      CQS2=AKHS02

!  END IF
  END SUBROUTINE SFCDIF3

  SUBROUTINE SFCDIF4(ILOC  ,JLOC  ,UX    ,VX     ,T1D  , &
                      P1D   ,PSFCPA,PBLH  ,DX     ,ZNT  , &
                      TSK   ,QX    ,ZLVL  ,IZ0TLND,QSFC , &
                      HFX   ,QFX   ,CM    ,CHS    ,CHS2 , &
                      CQS2  ,RMOL  ,UST   ,U10    ,V10)
                                                                                                                          
   USE MODULE_SF_SFCLAY                                                                                                   
   USE MODULE_MODEL_CONSTANTS                                                                                             
!-------------------------------------------------------------------                                                      
!  Compute surface drag coefficients CM for momentum and CH for heat                                                      
!  Joakim Refslund, 2011. Modified from YSU SFCLAY.                                                                                     
!-------------------------------------------------------------------                                                      
   IMPLICIT NONE                                                                                                          
!-------------------------------------------------------------------                                                      
! parameters                                                                                                              
   REAL,   PARAMETER     :: XKA=2.4E-5                                                                                    
   REAL,   PARAMETER     :: PRT=1.     !prandtl number                                                                    
                                                                                                                          
! input                                                                                                                   
   INTEGER,INTENT(IN )   :: ILOC                                                                                          
   INTEGER,INTENT(IN )   :: JLOC                                                                                          
                                                                                                                          
   REAL,   INTENT(IN )   :: PBLH      ! planetary boundary layer height                                                   
   REAL,   INTENT(IN )   :: TSK       ! skin temperature                                                                  
   REAL,   INTENT(IN )   :: PSFCPA    ! pressure in pascal                                                                
   REAL,   INTENT(IN )   :: P1D       !lowest model layer pressure (Pa)                                                      
   REAL,   INTENT(IN )   :: T1D       !lowest model layer temperature
!   REAL,   INTENT(IN )   :: QX        !water vapor mixing ratio (kg/kg)
   REAL,   INTENT(IN )   :: QX        !water vapor specific humidity (kg/kg)
   REAL,   INTENT(IN )   :: ZLVL      ! thickness of lowest full level layer
   REAL,   INTENT(IN )   :: HFX       ! sensible heat flux
   REAL,   INTENT(IN )   :: QFX       ! moisture flux
   REAL,   INTENT(IN )   :: DX        ! horisontal grid spacing
   REAL,   INTENT(IN )   :: UX
   REAL,   INTENT(IN )   :: VX
   REAL,   INTENT(IN )   :: ZNT                                                                                           
   REAL,   INTENT(INOUT ) :: QSFC

   REAL,   INTENT(INOUT) :: RMOL                                                                                          
   REAL,   INTENT(INOUT) :: UST                                                                                           
   REAL,   INTENT(INOUT) :: CHS2                                                                                          
   REAL,   INTENT(INOUT) :: CQS2                                                                                          
   REAL,   INTENT(INOUT) :: CHS                                                                                           
   REAL,   INTENT(INOUT) :: CM                
                                                                                                                          
! diagnostics out                                                                                                         
   REAL,   INTENT(OUT)   :: U10                                                                                           
   REAL,   INTENT(OUT)   :: V10                                                                                           
!   REAL,   INTENT(OUT)   :: TH2                                                                                           
!   REAL,   INTENT(OUT)   :: T2                                                                                            
!   REAL,   INTENT(OUT)   :: Q2                                                                                            
!   REAL,   INTENT(OUT)   :: QSFC                                                                                          
                                                                                                                          
! optional vars                                                                                                           
   INTEGER,OPTIONAL,INTENT(IN ) :: IZ0TLND                                                                                
                                                                                                                          
! local                                                                                                                   
   INTEGER :: REGIME  ! Stability regime
   REAL    :: ZA      ! Height of full-sigma level                                                                        
   REAL    :: THVX    ! Virtual potential temperature                                                                     
   REAL    :: ZQKL    ! Height of upper half level                                                                        
   REAL    :: ZQKLP1  ! Height of lower half level (surface)                                                              
   REAL    :: THX     ! Potential temperature                                                                             
   REAL    :: PSIH    ! similarity function for heat                                                                      
   REAL    :: PSIH2   ! Similarity function for heat 2m                                                                   
   REAL    :: PSIH10  ! Similarity function for heat 10m                                                                  
   REAL    :: PSIM    ! similarity function for momentum                                                                  
   REAL    :: PSIM2   ! Similarity function for momentum 2m                                                               
   REAL    :: PSIM10  ! Similarity function for momentum 10m                                                              
   REAL    :: DENOMQ  ! Denominator used for flux calc.                                                                   
   REAL    :: DENOMQ2 ! Denominator used for flux calc.                                                                   
   REAL    :: DENOMT2 ! Denominator used for flux calc.                                                                   
   REAL    :: WSPDI   ! Initial wind speed                                                                                
   REAL    :: GZ1OZ0  ! log(za/z0)                                                                                        
   REAL    :: GZ2OZ0  ! log(z2/z0)                                                                                        
   REAL    :: GZ10OZ0 ! log(z10/z0)                                                                                       
   REAL    :: RHOX    ! density                                                                                           
   REAL    :: GOVRTH  ! g/theta for stability L                                                                           
   REAL    :: TGDSA   ! tsk                                                                                               
!   REAL    :: SCR3    ! temporal variable -> input variable T1D                                                                                
   REAL    :: TVIR    ! temporal variable SRC4 -> TVIR                                                                                
   REAL    :: THGB    ! Potential temperature ground                                                                      
   REAL    :: PSFC    ! Surface pressure                                                                                  
   REAL    :: BR      ! bulk richardson number                                                                            
   REAL    :: CPM                                                                                           
   REAL    :: MOL                                                                                           
   REAL    :: ZOL  
   REAL    :: QGH    
   REAL    :: WSPD                                                                                          
                                                                                                                          
   INTEGER :: N,I,K,KK,L,NZOL,NK,NZOL2,NZOL10                                                                             
                                                                                                                          
   REAL    ::  PL,THCON,TVCON,E1                                                                                          
   REAL    ::  ZL,TSKV,DTHVDZ,DTHVM,VCONV,RZOL,RZOL2,RZOL10,ZOL2,ZOL10                                                    
   REAL    ::  DTG,PSIX,DTTHX,PSIX10,PSIT,PSIT2,PSIQ,PSIQ2,PSIQ10                                                         
   REAL    ::  FLUXC,VSGD,Z0Q,VISC,RESTAR,CZIL,RESTAR2                                                                    
!-------------------------------------------------------------------                                                      
   
   MOL = 1./RMOL
   ZL=0.01                                                                                                                
   PSFC=PSFCPA/1000.                                                                                                      
                                                                                                                          
! convert (tah or tgb = tsk) temperature to potential temperature.                                                                    
   TGDSA = TSK                                                                                                            
   THGB  = TSK*(P1000mb/PSFCPA)**RCP                                                                                      
                                                                                                                          
! store virtual, virtual potential and potential temperature
   PL    = P1D/1000.                                                                                                      
   THX   = T1D*(P1000mb*0.001/PL)**RCP                                                                                        
   THVX  = THX*(1.+EP_1*QX)                                                                                                   
   TVIR  = T1D*(1.+EP_1*QX)

! for land points QSFC can come from previous time step                                                                   
   !QSFC=EP_2*E1/(PSFC-E1)                                                                                
   IF (QSFC.LE.0.0) THEN
      !testing this
      E1=SVP1*EXP(SVP2*(TGDSA-SVPT0)/(TGDSA-SVP3))                                                                           
      QSFC=EP_2*E1/(PSFC-E1)                                                                                
      write(*,*) "JREF: IN SFCDIF4, QSFC WAS NEG. NOW = ",QSFC
   ENDIF

! qgh changed to use lowest-level air temp consistent with myjsfc change                                                  
! q2sat = qgh in lsm                                                                                                      
!jref: canres and esat is calculated in the loop so should that be changed??
!   QGH=EP_2*E1/(PL-E1)                                                                                                    
   CPM=CP*(1.+0.8*QX)                                                                                                     
                                                                                                                          
! compute the height of half-sigma levels above ground level                                                         
   !ZA=0.5*DZ8W                                                                                                   
   ZA = ZLVL
                                                                                                                          
! compute density and part of monin-obukhov length L                                                                      
   RHOX=PSFC*1000./(R_D*TVIR)                                                                                             
   GOVRTH=G/THX                                                                                                           
                                                                                                                          
! calculate bulk richardson no. of surface layer,                                                                         
! according to akb(1976), eq(12).                                                                                         
   GZ1OZ0=ALOG(ZA/ZNT)                                                                                                    
   GZ2OZ0=ALOG(2./ZNT)                                                                                                    
   GZ10OZ0=ALOG(10./ZNT)                                                                                                  
   WSPD=SQRT(UX*UX+VX*VX)                                                                                                 
                                                                                                                          
! virtual pot. temperature difference between input layer and lowest model layer                                              
   TSKV=THGB*(1.+EP_1*QSFC)                                                                                               
   DTHVDZ=(THVX-TSKV)                                                                                                     
                                                                                                                          
! convective velocity scale Vc and subgrid-scale velocity Vsg                                                             
! following Beljaars (1995, QJRMS) and Mahrt and Sun (1995, MWR)                                                          
!                                ... HONG Aug. 2001                                                                       
!                                                                                                                         
! VCONV = 0.25*sqrt(g/tskv*pblh(i)*dthvm)                                                                                 
! use Beljaars over land, old MM5 (Wyngaard) formula over water                                                           
                                                                                                                          
!jref:start commented out to see if stability is affected.                                                                                                                          
  FLUXC = MAX(HFX/RHOX/CP + EP_1*TSKV*QFX/RHOX,0.)                                                                                          
  VCONV = VCONVC*(G/TGDSA*PBLH*FLUXC)**.33                                                                                
!  VCONV = 0
!jref:end
                                                                                                                          
! Mahrt and Sun low-res correction                                                                                        
   VSGD = 0.32 * (max(dx/5000.-1.,0.))**.33                                                                               
   WSPD=SQRT(WSPD*WSPD+VCONV*VCONV+VSGD*VSGD)                                                                             
   WSPD=AMAX1(WSPD,0.1)                                                                                                   
   BR=GOVRTH*ZA*DTHVDZ/(WSPD*WSPD)                                                                                        
!  if previously unstable, do not let into regimes 1 and 2                                                                
   IF(MOL.LT.0.) BR=AMIN1(BR,0.0)                                                                                         
   RMOL=-GOVRTH*DTHVDZ*ZA*KARMAN                                                                                          
                                                                                                                          
!-----------------------------------------------------------------------                                                  
!     diagnose basic parameters for the appropriated stability class:                                                     
!                                                                                                                         
!     the stability classes are determined by br (bulk richardson no.)                                                    
!     and hol (height of pbl/monin-obukhov length).                                                                       
!                                                                                                                         
!     criteria for the classes are as follows:                                                                            
!                                                                                                                         
!        1. br .ge. 0.2;                                                                                                  
!               represents nighttime stable conditions (regime=1),                                                        
!                                                                                                                         
!        2. br .lt. 0.2 .and. br .gt. 0.0;                                                                                
!               represents damped mechanical turbulent conditions                                                         
!               (regime=2),                                                                                               
!                                                                                                                         
!        3. br .eq. 0.0                                                                                                   
!               represents forced convection conditions (regime=3),                                                       
!                                                                                                                         
!        4. br .lt. 0.0                                                                                                   
!               represents free convection conditions (regime=4).                                                         
!                                                                                                                         
!-----------------------------------------------------------------------                                                  
                                                                                                                          
   IF (BR.GE.0.2) REGIME=1                                                                                                
   IF (BR.LT.0.2 .AND. BR.GT.0.0) REGIME=2                                                                                
   IF (BR.EQ.0.0) REGIME=3                                                                                                
   IF (BR.LT.0.0) REGIME=4                                                                                                
                                                                                                                          
   SELECT CASE(REGIME)                                                                                                    
     CASE(1) ! class 1; stable (nighttime) conditions:                                                                    
       PSIM=-10.*GZ1OZ0                                                                                                   
! lower limit on psi in stable conditions                                                                                 
       PSIM=AMAX1(PSIM,-10.)                                                                                              
       PSIH=PSIM                                                                                                          
       PSIM10=10./ZA*PSIM                                                                                                 
       PSIM10=AMAX1(PSIM10,-10.)                                                                                          
       PSIH10=PSIM10                                                                                                      
       PSIM2=2./ZA*PSIM                                                                                                   
       PSIM2=AMAX1(PSIM2,-10.)                                                                                            
       PSIH2=PSIM2                                                                                                        
                                                                                                                          
! 1.0 over Monin-Obukhov length                                                                                           
       IF(UST.LT.0.01)THEN                                                                                                
          RMOL=BR*GZ1OZ0 !ZA/L                                                                                            
       ELSE                                                                                                               
          RMOL=KARMAN*GOVRTH*ZA*MOL/(UST*UST) !ZA/L                                                                       
       ENDIF                                                                                                              
       RMOL=AMIN1(RMOL,9.999) ! ZA/L                                                                                      
       RMOL = RMOL/ZA !1.0/L                                                                                              
                                                                                                                          
     CASE(2) ! class 2; damped mechanical turbulence:                                                                     
       PSIM=-5.0*BR*GZ1OZ0/(1.1-5.0*BR)                                                                                   
! lower limit on psi in stable conditions                                                                                 
       PSIM=AMAX1(PSIM,-10.)                                                                                              
! AKB(1976), EQ(16).                                                                                                      
       PSIH=PSIM                                                                                                          
       PSIM10=10./ZA*PSIM                                                                                                 
       PSIM10=AMAX1(PSIM10,-10.)                                                                                          
       PSIH10=PSIM10                                                                                                      
       PSIM2=2./ZA*PSIM                                                                                                   
       PSIM2=AMAX1(PSIM2,-10.)                                                                                            
       PSIH2=PSIM2                                                                                                        
                                                                                                                          
     ! Linear form: PSIM = -0.5*ZA/L; e.g, see eqn 16 of                                                                  
     ! Blackadar, Modeling the nocturnal boundary layer, Preprints,                                                       
     ! Third Symposium on Atmospheric Turbulence Diffusion and Air Quality,                                               
     ! Raleigh, NC, 1976                                                                                                  
       ZOL = BR*GZ1OZ0/(1.00001-5.0*BR)                                                                                   
                                                                                                                          
       IF ( ZOL .GT. 0.5 ) THEN ! linear form ok                                                                          
        ! Holtslag and de Bruin, J. App. Meteor 27, 689-704, 1988;                                                        
        ! see also, Launiainen, Boundary-Layer Meteor 76,165-179, 1995                                                    
        ! Eqn (8) of Launiainen, 1995                                                                                     
          ZOL = ( 1.89*GZ1OZ0 + 44.2 ) * BR*BR    &                                                                       
               + ( 1.18*GZ1OZ0 - 1.37 ) * BR                                                                              
          ZOL=AMIN1(ZOL,9.999)                                                                                            
       END IF                                                                                                             
                                                                                                                          
! 1.0 over Monin-Obukhov length                                                                                           
       RMOL= ZOL/ZA                                                                                                       
                                                                                                                          
     CASE(3)  ! class 3; forced convection:                                                                               
       PSIM=0.0                                                                                                           
       PSIH=PSIM                                                                                                          
       PSIM10=0.                                                                                                          
       PSIH10=PSIM10                                                                                                      
       PSIM2=0.                                                                                                           
       PSIH2=PSIM2                                                                                                        
       IF(UST.LT.0.01)THEN                                                                                                
         ZOL=BR*GZ1OZ0                                                                                                    
       ELSE                                                                                                               
         ZOL=KARMAN*GOVRTH*ZA*MOL/(UST*UST)                                                                               
       ENDIF                                                                                                              
                                                                                                                          
       RMOL = ZOL/ZA                                                                                                      
                                                                                                                          
     CASE(4) ! class 4; free convection:                                                                                  
       IF(UST.LT.0.01)THEN                                                                                                
         ZOL=BR*GZ1OZ0                                                                                                    
       ELSE                                                                                                               
         ZOL=KARMAN*GOVRTH*ZA*MOL/(UST*UST)                                                                               
       ENDIF                                                                                                              
       ZOL10=10./ZA*ZOL                                                                                                   
       ZOL2=2./ZA*ZOL                                                                                                     
       ZOL=AMIN1(ZOL,0.)                                                                                                  
       ZOL=AMAX1(ZOL,-9.9999)                                                                                             
       ZOL10=AMIN1(ZOL10,0.)                                                                                              
       ZOL10=AMAX1(ZOL10,-9.9999)                                                                                         
       ZOL2=AMIN1(ZOL2,0.)                                                                                                
       ZOL2=AMAX1(ZOL2,-9.9999)                                                                                           
       NZOL=INT(-ZOL*100.)                                                                                                
       RZOL=-ZOL*100.-NZOL                                                                                                
       NZOL10=INT(-ZOL10*100.)                                                                                            
       RZOL10=-ZOL10*100.-NZOL10                                                                                          
       NZOL2=INT(-ZOL2*100.)                                                                                              
       RZOL2=-ZOL2*100.-NZOL2                                                                                             
       PSIM=PSIMTB(NZOL)+RZOL*(PSIMTB(NZOL+1)-PSIMTB(NZOL))                                                               
       PSIH=PSIHTB(NZOL)+RZOL*(PSIHTB(NZOL+1)-PSIHTB(NZOL))                                                               
       PSIM10=PSIMTB(NZOL10)+RZOL10*(PSIMTB(NZOL10+1)-PSIMTB(NZOL10))                                                     
       PSIH10=PSIHTB(NZOL10)+RZOL10*(PSIHTB(NZOL10+1)-PSIHTB(NZOL10))                                                     
       PSIM2=PSIMTB(NZOL2)+RZOL2*(PSIMTB(NZOL2+1)-PSIMTB(NZOL2))                                                          
       PSIH2=PSIHTB(NZOL2)+RZOL2*(PSIHTB(NZOL2+1)-PSIHTB(NZOL2))                                                          
                                                                                                                          
! limit psih and psim in the case of thin layers and high roughness                                                       
! this prevents denominator in fluxes from getting too small                                                              
!       PSIH=AMIN1(PSIH,0.9*GZ1OZ0)                                                                                       
!       PSIM=AMIN1(PSIM,0.9*GZ1OZ0)                                                                                       
       PSIH=AMIN1(PSIH,0.9*GZ1OZ0)                                                                                        
       PSIM=AMIN1(PSIM,0.9*GZ1OZ0)                                                                                        
       PSIH2=AMIN1(PSIH2,0.9*GZ2OZ0)                                                                                      
       PSIM10=AMIN1(PSIM10,0.9*GZ10OZ0)                                                                                   
! AHW: mods to compute ck, cd                                                                                             
       PSIH10=AMIN1(PSIH10,0.9*GZ10OZ0)                                                                                   
                                                                                                                          
       RMOL = ZOL/ZA                                                                                                      
                                                                                                                          
   END SELECT ! stability regime done                                                                                     
                                                                                                                          
! compute the frictional velocity: ZA(1982) EQS(2.60),(2.61).                                                             
   DTG=THX-THGB                                                                                                           
   PSIX=GZ1OZ0-PSIM                                                                                                       
   PSIX10=GZ10OZ0-PSIM10                                                                                                  
                                                                                                                          
! lower limit added to prevent large flhc in soil model                                                                   
! activates in unstable conditions with thin layers or high z0                                                            
   PSIT=AMAX1(GZ1OZ0-PSIH,2.) !does this still apply???? jref                                                             
   PSIQ=ALOG(KARMAN*UST*ZA/XKA+ZA/ZL)-PSIH                                                                                
   PSIT2=GZ2OZ0-PSIH2                                                                                                     
   PSIQ2=ALOG(KARMAN*UST*2./XKA+2./ZL)-PSIH2                                                                              
! AHW: mods to compute ck, cd                                                                                             
   PSIQ10=ALOG(KARMAN*UST*10./XKA+10./ZL)-PSIH10                                                                          

!jref:start - commented out since these values can be produced by sfclay routine   
!   IF(PRESENT(ck) .and. PRESENT(cd) .and. PRESENT(cka) .and. PRESENT(cda)) THEN                                           
!      Ck=(karman/psix10)*(karman/psiq10)                                                                                  
!      Cd=(karman/psix10)*(karman/psix10)                                                                                  
!      Cka=(karman/psix)*(karman/psiq)                                                                                     
!      Cda=(karman/psix)*(karman/psix)                                                                                     
!   ENDIF                                                                                                                  

!   WRITE(*,*) "KARMAN=",KARMAN
!   WRITE(*,*) "UST=",UST
!   WRITE(*,*) "XKA=",XKA
!   WRITE(*,*) "ZA =",ZA
!   WRITE(*,*) "ZL =",ZL
!   WRITE(*,*) "PSIH=",PSIH
!   WRITE(*,*) "PSIQ=",PSIQ,"PSIT=",PSIT

   IF ( PRESENT(IZ0TLND) ) THEN                                                                                           
      IF ( IZ0TLND.EQ.1 ) THEN                                                                                            
         ZL=ZNT                                                                                                           
!       czil related changes for land                                                                                     
         VISC=(1.32+0.009*(T1D-273.15))*1.E-5                                                                            
         RESTAR=UST*ZL/VISC                                                                                               
!       modify CZIL according to Chen & Zhang, 2009                                                                       
                                                                                                                          
         CZIL = 10.0 ** ( -0.40 * ( ZL / 0.07 ) )                                                                         
                                                                                                                          
         PSIT=GZ1OZ0-PSIH+CZIL*KARMAN*SQRT(RESTAR)                                                                        
         PSIQ=GZ1OZ0-PSIH+CZIL*KARMAN*SQRT(RESTAR)                                                                        
         PSIT2=GZ2OZ0-PSIH2+CZIL*KARMAN*SQRT(RESTAR)                                                                      
         PSIQ2=GZ2OZ0-PSIH2+CZIL*KARMAN*SQRT(RESTAR)                                                                      
      ENDIF                                                                                                               
   ENDIF                                                                                                                  
                                                                                                                          

! to prevent oscillations average with old value                                                                          
   UST=0.5*UST+0.5*KARMAN*WSPD/PSIX                                                                                       
   UST=AMAX1(UST,0.1)                                                                                                     
!jref: should this be converted to RMOL???
   MOL=KARMAN*DTG/PSIT/PRT                                                                                                
   DENOMQ=PSIQ                                                                                                            
   DENOMQ2=PSIQ2                                                                                                          
   DENOMT2=PSIT2                                                                                                          
!   WRITE(*,*) "ILOC,JLOC=",ILOC,JLOC,"DENOMQ=",DENOMQ
!   WRITE(*,*) "UST=",UST,"PSIT=",PSIT
!   call wrf_error_fatal("stop in sfcdif4")
                                                                                                                          
! calculate exchange coefficients                                                                                         
!jref: start exchange coefficient for momentum
   CM =KARMAN*KARMAN/(PSIX*PSIX)
!jref:end
   CHS=UST*KARMAN/DENOMQ                                                                                                  
!        GZ2OZ0=ALOG(2./ZNT)                                                                                              
!        PSIM2=-10.*GZ2OZ0                                                                                                
!        PSIM2=AMAX1(PSIM2,-10.)                                                                                          
!        PSIH2=PSIM2                                                                                                      
   CQS2=UST*KARMAN/DENOMQ2                                                                                                
   CHS2=UST*KARMAN/DENOMT2                                                                                                
! jref: in last iteration calculate diagnostics                                                                           
                                                                                                                         
   U10=UX*PSIX10/PSIX                                                                                                     
   V10=VX*PSIX10/PSIX                                                                                                     
                                                                                                                          
! jref: check the following for correct calculation                                                                       
!   TH2=THGB+DTG*PSIT2/PSIT                                                                                               
!   Q2=QSFC+(QX-QSFC)*PSIQ2/PSIQ                                                                                          
!   T2 = TH2*(PSFCPA/P1000mb)**RCP                                                                                        
                                                                                                                          
   END SUBROUTINE SFCDIF4   
