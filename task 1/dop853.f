      SUBROUTINE DOP853(N,FCN,X,Y,XEND,
     &                  RTOL,ATOL,ITOL,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C ----------------------------------------------------------
C     NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER
C     ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y).
C     THIS IS AN EXPLICIT RUNGE-KUTTA METHOD OF ORDER 8(5,3)  
C     DUE TO DORMAND & PRINCE (WITH STEPSIZE CONTROL AND
C     DENSE OUTPUT)
C
C     AUTHORS: E. HAIRER AND G. WANNER
C              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
C              CH-1211 GENEVE 24, SWITZERLAND 
C              E-MAIL:  Ernst.Hairer@math.unige.ch
C                       Gerhard.Wanner@math.unige.ch
C     
C     THIS CODE IS DESCRIBED IN:
C         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
C         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION. 
C         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS, 
C         SPRINGER-VERLAG (1993)
C      
C     VERSION OF APRIL 25, 1996
C         SMALL CORRECTIONS ON JULY 12, 2000
C
C     INPUT PARAMETERS  
C     ----------------  
C     N           DIMENSION OF THE SYSTEM 
C
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
C                 VALUE OF F(X,Y):
C                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(N),F(N)
C                    F(1)=...   ETC.
C
C     X           INITIAL X-VALUE
C
C     Y(N)        INITIAL VALUES FOR Y
C
C     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
C
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
C                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
C                 ATOL SHOULD BE STRICTLY POSITIVE (POSSIBLY VERY SMALL)
C
C     ITOL        SWITCH FOR RTOL AND ATOL:
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
C                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
C                     RTOL(I)*ABS(Y(I))+ATOL(I).
C
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION. 
C                 IF IOUT.GE.1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
C                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. 
C                 IT MUST HAVE THE FORM
C                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,
C                                       RPAR,IPAR,IRTRN)
C                    DIMENSION Y(N),CON(8*ND),ICOMP(ND)
C                    ....  
C                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
C                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
C                    THE FIRST GRID-POINT).
C                 "XOLD" IS THE PRECEEDING GRID-POINT.
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, DOP853 WILL RETURN TO THE CALLING PROGRAM.
C                    IF THE NUMERICAL SOLUTION IS ALTERED IN SOLOUT,
C                    SET  IRTRN = 2
C           
C          -----  CONTINUOUS OUTPUT: -----
C                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
C                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
C                 THE FUNCTION
C                        >>>   CONTD8(I,S,CON,ICOMP,ND)   <<<
C                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
C                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
C                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
C           
C     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
C                    IOUT=0: SUBROUTINE IS NEVER CALLED
C                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
C                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT
C                            (IN THIS CASE WORK(5) MUST BE SPECIFIED)
C
C     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
C                 WORK(1),...,WORK(20) SERVE AS PARAMETERS FOR THE CODE.
C                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING.
C                 "LWORK" MUST BE AT LEAST  11*N+8*NRDENS+20
C                 WHERE  NRDENS = IWORK(5)
C
C     LWORK       DECLARED LENGHT OF ARRAY "WORK".
C
C     IWORK       INTEGER WORKING SPACE OF LENGHT "LIWORK".
C                 IWORK(1),...,IWORK(20) SERVE AS PARAMETERS FOR THE CODE.
C                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING.
C                 "LIWORK" MUST BE AT LEAST NRDENS+20 .
C
C     LIWORK      DECLARED LENGHT OF ARRAY "IWORK".
C
C     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH  
C                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
C                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES. 
C
C-----------------------------------------------------------------------
C 
C     SOPHISTICATED SETTING OF PARAMETERS
C     -----------------------------------
C              SEVERAL PARAMETERS (WORK(1),...,IWORK(1),...) ALLOW
C              TO ADAPT THE CODE TO THE PROBLEM AND TO THE NEEDS OF
C              THE USER. FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES.
C
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 2.3D-16.
C
C    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
C              DEFAULT 0.9D0.
C
C    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION
C              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C                 WORK(3) <= HNEW/HOLD <= WORK(4)
C              DEFAULT VALUES: WORK(3)=0.333D0, WORK(4)=6.D0
C
C    WORK(5)   IS THE "BETA" FOR STABILIZED STEP SIZE CONTROL
C              (SEE SECTION IV.2). POSITIVE VALUES OF BETA ( <= 0.04 )
C              MAKE THE STEP SIZE CONTROL MORE STABLE.
C              NEGATIVE WORK(5) PROVOKE BETA=0.
C              DEFAULT 0.0D0.
C
C    WORK(6)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
C
C    WORK(7)   INITIAL STEP SIZE, FOR WORK(7)=0.D0 AN INITIAL GUESS
C              IS COMPUTED WITH HELP OF THE FUNCTION HINIT
C
C    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000.
C
C    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS
C              IF IWORK(2).EQ.1  METHOD DOP853 OF DORMAND AND PRINCE
C              (SECTION II.6).
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=1.
C
C    IWORK(3)  SWITCH FOR PRINTING ERROR MESSAGES
C              IF IWORK(3).LT.0 NO MESSAGES ARE BEING PRINTED
C              IF IWORK(3).GT.0 MESSAGES ARE PRINTED WITH
C              WRITE (IWORK(3),*) ...  
C              DEFAULT VALUE (FOR IWORK(3)=0) IS IWORK(3)=6
C
C    IWORK(4)  TEST FOR STIFFNESS IS ACTIVATED AFTER STEP NUMBER
C              J*IWORK(4) (J INTEGER), PROVIDED IWORK(4).GT.0.
C              FOR NEGATIVE IWORK(4) THE STIFFNESS TEST IS
C              NEVER ACTIVATED; DEFAULT VALUE IS IWORK(4)=1000
C
C    IWORK(5)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT
C              IS REQUIRED; DEFAULT VALUE IS IWORK(5)=0;
C              FOR   0 < NRDENS < N   THE COMPONENTS (FOR WHICH DENSE
C              OUTPUT IS REQUIRED) HAVE TO BE SPECIFIED IN
C              IWORK(21),...,IWORK(NRDENS+20);
C              FOR  NRDENS=N  THIS IS DONE BY THE CODE.
C
C----------------------------------------------------------------------
C
C     OUTPUT PARAMETERS 
C     ----------------- 
C     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
C                 (AFTER SUCCESSFUL RETURN X=XEND).
C
C     Y(N)        NUMERICAL SOLUTION AT X
C 
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
C
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
C                   IDID= 1  COMPUTATION SUCCESSFUL,
C                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
C                   IDID=-1  INPUT IS NOT CONSISTENT,
C                   IDID=-2  LARGER NMAX IS NEEDED,
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL.
C                   IDID=-4  PROBLEM IS PROBABLY STIFF (INTERRUPTED).
C
C   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS
C   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS
C   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS
C   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
C                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
C-----------------------------------------------------------------------
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C          DECLARATIONS 
C *** *** *** *** *** *** *** *** *** *** *** *** ***
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),ATOL(*),RTOL(*),WORK(LWORK),IWORK(LIWORK)
      DIMENSION RPAR(*),IPAR(*)
      LOGICAL ARRET
      EXTERNAL FCN,SOLOUT
C *** *** *** *** *** *** ***
C        SETTING THE PARAMETERS 
C *** *** *** *** *** *** ***
      NFCN=0
      NSTEP=0
      NACCPT=0
      NREJCT=0
      ARRET=.FALSE.
C -------- IPRINT FOR MONITORING THE PRINTING
      IF(IWORK(3).EQ.0)THEN
         IPRINT=6
      ELSE
         IPRINT=IWORK(3)
      END IF
C -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- 
      IF(IWORK(1).EQ.0)THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(1)
         IF(NMAX.LE.0)THEN
            IF (IPRINT.GT.0) WRITE(IPRINT,*)
     &          ' WRONG INPUT IWORK(1)=',IWORK(1)
            ARRET=.TRUE.
         END IF
      END IF
C -------- METH   COEFFICIENTS OF THE METHOD
      IF(IWORK(2).EQ.0)THEN
         METH=1
      ELSE
         METH=IWORK(2)
         IF(METH.LE.0.OR.METH.GE.4)THEN
            IF (IPRINT.GT.0) WRITE(IPRINT,*)
     &          ' CURIOUS INPUT IWORK(2)=',IWORK(2)
            ARRET=.TRUE.
         END IF
      END IF  
C -------- NSTIFF   PARAMETER FOR STIFFNESS DETECTION  
      NSTIFF=IWORK(4) 
      IF (NSTIFF.EQ.0) NSTIFF=1000
      IF (NSTIFF.LT.0) NSTIFF=NMAX+10
C -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS
      NRDENS=IWORK(5)
      IF(NRDENS.LT.0.OR.NRDENS.GT.N)THEN
         IF (IPRINT.GT.0) WRITE(IPRINT,*)
     &           ' CURIOUS INPUT IWORK(5)=',IWORK(5)
         ARRET=.TRUE.
      ELSE
         IF(NRDENS.GT.0.AND.IOUT.LT.2)THEN
            IF (IPRINT.GT.0) WRITE(IPRINT,*)
     &       ' WARNING: PUT IOUT=2 FOR DENSE OUTPUT '
         END IF 
         IF (NRDENS.EQ.N) THEN
            DO I=1,NRDENS
               IWORK(I+20)=I
            END DO
         END IF
      END IF       
C -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0  
      IF(WORK(1).EQ.0.D0)THEN
         UROUND=2.3D-16
      ELSE
         UROUND=WORK(1)
         IF(UROUND.LE.1.D-35.OR.UROUND.GE.1.D0)THEN
            IF (IPRINT.GT.0) WRITE(IPRINT,*)
     &        ' WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:',WORK(1)
            ARRET=.TRUE.
         END IF
      END IF
C -------  SAFETY FACTOR -------------
      IF(WORK(2).EQ.0.D0)THEN
         SAFE=0.9D0
      ELSE
         SAFE=WORK(2)
         IF(SAFE.GE.1.D0.OR.SAFE.LE.1.D-4)THEN
            IF (IPRINT.GT.0) WRITE(IPRINT,*)
     &          ' CURIOUS INPUT FOR SAFETY FACTOR WORK(2)=',WORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
      IF(WORK(3).EQ.0.D0)THEN
         FAC1=0.333D0
      ELSE
         FAC1=WORK(3)
      END IF
      IF(WORK(4).EQ.0.D0)THEN
         FAC2=6.D0
      ELSE
         FAC2=WORK(4)
      END IF
C --------- BETA FOR STEP CONTROL STABILIZATION -----------
      IF(WORK(5).EQ.0.D0)THEN
         BETA=0.0D0
      ELSE
         IF(WORK(5).LT.0.D0)THEN
            BETA=0.D0
         ELSE
            BETA=WORK(5)
            IF(BETA.GT.0.2D0)THEN
               IF (IPRINT.GT.0) WRITE(IPRINT,*)
     &          ' CURIOUS INPUT FOR BETA: WORK(5)=',WORK(5)
            ARRET=.TRUE.
         END IF
         END IF
      END IF
C -------- MAXIMAL STEP SIZE
      IF(WORK(6).EQ.0.D0)THEN
         HMAX=XEND-X
      ELSE
         HMAX=WORK(6)
      END IF
C -------- INITIAL STEP SIZE
      H=WORK(7)
C ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
      IEK1=21
      IEK2=IEK1+N
      IEK3=IEK2+N
      IEK4=IEK3+N
      IEK5=IEK4+N
      IEK6=IEK5+N
      IEK7=IEK6+N
      IEK8=IEK7+N
      IEK9=IEK8+N
      IEK10=IEK9+N
      IEY1=IEK10+N
      IECO=IEY1+N
C ------ TOTAL STORAGE REQUIREMENT -----------
      ISTORE=IECO+8*NRDENS-1
      IF(ISTORE.GT.LWORK)THEN
        IF (IPRINT.GT.0) WRITE(IPRINT,*)
     &   ' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
        ARRET=.TRUE.
      END IF
      ICOMP=21
      ISTORE=ICOMP+NRDENS-1
      IF(ISTORE.GT.LIWORK)THEN
        IF (IPRINT.GT.0) WRITE(IPRINT,*)
     &   ' INSUFFICIENT STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
        ARRET=.TRUE.
      END IF
C -------- WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF
C -------- CALL TO CORE INTEGRATOR ------------
      CALL DP86CO(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,IPRINT,
     &   SOLOUT,IOUT,IDID,NMAX,UROUND,METH,NSTIFF,SAFE,BETA,FAC1,FAC2,
     &   WORK(IEK1),WORK(IEK2),WORK(IEK3),WORK(IEK4),WORK(IEK5),
     &   WORK(IEK6),WORK(IEK7),WORK(IEK8),WORK(IEK9),WORK(IEK10),
     &   WORK(IEY1),WORK(IECO),IWORK(ICOMP),NRDENS,RPAR,IPAR,
     &   NFCN,NSTEP,NACCPT,NREJCT)
      WORK(7)=H
      IWORK(17)=NFCN
      IWORK(18)=NSTEP
      IWORK(19)=NACCPT
      IWORK(20)=NREJCT
C ----------- RETURN -----------
      RETURN
      END
C
C
C
C  ----- ... AND HERE IS THE CORE INTEGRATOR  ----------
C
      SUBROUTINE DP86CO(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,IPRINT,
     &   SOLOUT,IOUT,IDID,NMAX,UROUND,METH,NSTIFF,SAFE,BETA,FAC1,FAC2,
     &   K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,Y1,CONT,ICOMP,NRD,RPAR,IPAR,
     &   NFCN,NSTEP,NACCPT,NREJCT)
C ----------------------------------------------------------
C     CORE INTEGRATOR FOR DOP853
C     PARAMETERS SAME AS IN DOP853 WITH WORKSPACE ADDED 
C ---------------------------------------------------------- 
C         DECLARATIONS 
C ---------------------------------------------------------- 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter (
     &  c2  = 0.526001519587677318785587544488D-01,
     &  c3  = 0.789002279381515978178381316732D-01,
     &  c4  = 0.118350341907227396726757197510D+00,
     &  c5  = 0.281649658092772603273242802490D+00,
     &  c6  = 0.333333333333333333333333333333D+00,
     &  c7  = 0.25D+00,
     &  c8  = 0.307692307692307692307692307692D+00,
     &  c9  = 0.651282051282051282051282051282D+00,
     &  c10 = 0.6D+00,
     &  c11 = 0.857142857142857142857142857142D+00,
     &  c14 = 0.1D+00,
     &  c15 = 0.2D+00,
     &  c16 = 0.777777777777777777777777777778D+00)
      parameter (
     &  b1 =   5.42937341165687622380535766363D-2,
     &  b6 =   4.45031289275240888144113950566D0,
     &  b7 =   1.89151789931450038304281599044D0,
     &  b8 =  -5.8012039600105847814672114227D0,
     &  b9 =   3.1116436695781989440891606237D-1,
     &  b10 = -1.52160949662516078556178806805D-1,
     &  b11 =  2.01365400804030348374776537501D-1,
     &  b12 =  4.47106157277725905176885569043D-2)
      parameter (
     &  bhh1 = 0.244094488188976377952755905512D+00,
     &  bhh2 = 0.733846688281611857341361741547D+00,
     &  bhh3 = 0.220588235294117647058823529412D-01)
      parameter (
     &  er 1 =  0.1312004499419488073250102996D-01,
     &  er 6 = -0.1225156446376204440720569753D+01,
     &  er 7 = -0.4957589496572501915214079952D+00,
     &  er 8 =  0.1664377182454986536961530415D+01,
     &  er 9 = -0.3503288487499736816886487290D+00,
     &  er10 =  0.3341791187130174790297318841D+00,
     &  er11 =  0.8192320648511571246570742613D-01,
     &  er12 = -0.2235530786388629525884427845D-01)
      parameter (
     &  a21 =    5.26001519587677318785587544488D-2,
     &  a31 =    1.97250569845378994544595329183D-2,
     &  a32 =    5.91751709536136983633785987549D-2,
     &  a41 =    2.95875854768068491816892993775D-2,
     &  a43 =    8.87627564304205475450678981324D-2,
     &  a51 =    2.41365134159266685502369798665D-1,
     &  a53 =   -8.84549479328286085344864962717D-1,
     &  a54 =    9.24834003261792003115737966543D-1,
     &  a61 =    3.7037037037037037037037037037D-2,
     &  a64 =    1.70828608729473871279604482173D-1,
     &  a65 =    1.25467687566822425016691814123D-1,
     &  a71 =    3.7109375D-2,
     &  a74 =    1.70252211019544039314978060272D-1,
     &  a75 =    6.02165389804559606850219397283D-2,
     &  a76 =   -1.7578125D-2)
      parameter (
     &  a81 =    3.70920001185047927108779319836D-2,
     &  a84 =    1.70383925712239993810214054705D-1,
     &  a85 =    1.07262030446373284651809199168D-1,
     &  a86 =   -1.53194377486244017527936158236D-2,
     &  a87 =    8.27378916381402288758473766002D-3,
     &  a91 =    6.24110958716075717114429577812D-1,
     &  a94 =   -3.36089262944694129406857109825D0,
     &  a95 =   -8.68219346841726006818189891453D-1,
     &  a96 =    2.75920996994467083049415600797D1,
     &  a97 =    2.01540675504778934086186788979D1,
     &  a98 =   -4.34898841810699588477366255144D1,
     &  a101 =   4.77662536438264365890433908527D-1,
     &  a104 =  -2.48811461997166764192642586468D0,
     &  a105 =  -5.90290826836842996371446475743D-1,
     &  a106 =   2.12300514481811942347288949897D1,
     &  a107 =   1.52792336328824235832596922938D1,
     &  a108 =  -3.32882109689848629194453265587D1,
     &  a109 =  -2.03312017085086261358222928593D-2)
      parameter (
     &  a111 =  -9.3714243008598732571704021658D-1,
     &  a114 =   5.18637242884406370830023853209D0,
     &  a115 =   1.09143734899672957818500254654D0,
     &  a116 =  -8.14978701074692612513997267357D0,
     &  a117 =  -1.85200656599969598641566180701D1,
     &  a118 =   2.27394870993505042818970056734D1,
     &  a119 =   2.49360555267965238987089396762D0,
     &  a1110 = -3.0467644718982195003823669022D0,
     &  a121 =   2.27331014751653820792359768449D0,
     &  a124 =  -1.05344954667372501984066689879D1,
     &  a125 =  -2.00087205822486249909675718444D0,
     &  a126 =  -1.79589318631187989172765950534D1,
     &  a127 =   2.79488845294199600508499808837D1,
     &  a128 =  -2.85899827713502369474065508674D0,
     &  a129 =  -8.87285693353062954433549289258D0,
     &  a1210 =  1.23605671757943030647266201528D1,
     &  a1211 =  6.43392746015763530355970484046D-1)
      parameter (
     &  a141 =  5.61675022830479523392909219681D-2,
     &  a147 =  2.53500210216624811088794765333D-1,
     &  a148 = -2.46239037470802489917441475441D-1,
     &  a149 = -1.24191423263816360469010140626D-1,
     &  a1410 =  1.5329179827876569731206322685D-1,
     &  a1411 =  8.20105229563468988491666602057D-3,
     &  a1412 =  7.56789766054569976138603589584D-3,
     &  a1413 = -8.298D-3)
      parameter (
     &  a151 =  3.18346481635021405060768473261D-2,
     &  a156 =  2.83009096723667755288322961402D-2,
     &  a157 =  5.35419883074385676223797384372D-2,
     &  a158 = -5.49237485713909884646569340306D-2,
     &  a1511 = -1.08347328697249322858509316994D-4,
     &  a1512 =  3.82571090835658412954920192323D-4,
     &  a1513 = -3.40465008687404560802977114492D-4,
     &  a1514 =  1.41312443674632500278074618366D-1,
     &  a161 = -4.28896301583791923408573538692D-1,
     &  a166 = -4.69762141536116384314449447206D0,
     &  a167 =  7.68342119606259904184240953878D0,
     &  a168 =  4.06898981839711007970213554331D0,
     &  a169 =  3.56727187455281109270669543021D-1,
     &  a1613 = -1.39902416515901462129418009734D-3,
     &  a1614 =  2.9475147891527723389556272149D0,
     &  a1615 = -9.15095847217987001081870187138D0)
      parameter (
     &  d41  = -0.84289382761090128651353491142D+01,
     &  d46  =  0.56671495351937776962531783590D+00,
     &  d47  = -0.30689499459498916912797304727D+01,
     &  d48  =  0.23846676565120698287728149680D+01,
     &  d49  =  0.21170345824450282767155149946D+01,
     &  d410 = -0.87139158377797299206789907490D+00,
     &  d411 =  0.22404374302607882758541771650D+01,
     &  d412 =  0.63157877876946881815570249290D+00,
     &  d413 = -0.88990336451333310820698117400D-01,
     &  d414 =  0.18148505520854727256656404962D+02,
     &  d415 = -0.91946323924783554000451984436D+01,
     &  d416 = -0.44360363875948939664310572000D+01)
      parameter (
     &  d51  =  0.10427508642579134603413151009D+02,
     &  d56  =  0.24228349177525818288430175319D+03,
     &  d57  =  0.16520045171727028198505394887D+03,
     &  d58  = -0.37454675472269020279518312152D+03,
     &  d59  = -0.22113666853125306036270938578D+02,
     &  d510 =  0.77334326684722638389603898808D+01,
     &  d511 = -0.30674084731089398182061213626D+02,
     &  d512 = -0.93321305264302278729567221706D+01,
     &  d513 =  0.15697238121770843886131091075D+02,
     &  d514 = -0.31139403219565177677282850411D+02,
     &  d515 = -0.93529243588444783865713862664D+01,
     &  d516 =  0.35816841486394083752465898540D+02)
      parameter (
     &  d61 =  0.19985053242002433820987653617D+02,
     &  d66 = -0.38703730874935176555105901742D+03,
     &  d67 = -0.18917813819516756882830838328D+03,
     &  d68 =  0.52780815920542364900561016686D+03,
     &  d69 = -0.11573902539959630126141871134D+02,
     &  d610 =  0.68812326946963000169666922661D+01,
     &  d611 = -0.10006050966910838403183860980D+01,
     &  d612 =  0.77771377980534432092869265740D+00,
     &  d613 = -0.27782057523535084065932004339D+01,
     &  d614 = -0.60196695231264120758267380846D+02,
     &  d615 =  0.84320405506677161018159903784D+02,
     &  d616 =  0.11992291136182789328035130030D+02)
      parameter (
     &  d71  = -0.25693933462703749003312586129D+02,
     &  d76  = -0.15418974869023643374053993627D+03,
     &  d77  = -0.23152937917604549567536039109D+03,
     &  d78  =  0.35763911791061412378285349910D+03,
     &  d79  =  0.93405324183624310003907691704D+02,
     &  d710 = -0.37458323136451633156875139351D+02,
     &  d711 =  0.10409964950896230045147246184D+03,
     &  d712 =  0.29840293426660503123344363579D+02,
     &  d713 = -0.43533456590011143754432175058D+02,
     &  d714 =  0.96324553959188282948394950600D+02,
     &  d715 = -0.39177261675615439165231486172D+02,
     &  d716 = -0.14972683625798562581422125276D+03)
      DOUBLE PRECISION Y(N),Y1(N),K1(N),K2(N),K3(N),K4(N),K5(N),K6(N)
      DOUBLE PRECISION K7(N),K8(N),K9(N),K10(N),ATOL(*),RTOL(*)     
      DIMENSION CONT(8*NRD),ICOMP(NRD),RPAR(*),IPAR(*)
      LOGICAL REJECT,LAST 
      EXTERNAL FCN
      COMMON /CONDO8/XOLD,HOUT
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** *** 
      FACOLD=1.D-4  
      EXPO1=1.d0/8.d0-BETA*0.2D0
      FACC1=1.D0/FAC1
      FACC2=1.D0/FAC2
      POSNEG=SIGN(1.D0,XEND-X) 
C --- INITIAL PREPARATIONS   
      ATOLI=ATOL(1)
      RTOLI=RTOL(1)    
      LAST=.FALSE. 
      HLAMB=0.D0
      IASTI=0
      CALL FCN(N,X,Y,K1,RPAR,IPAR)
      HMAX=ABS(HMAX)     
      IORD=8  
      IF (H.EQ.0.D0) H=HINIT(N,FCN,X,Y,XEND,POSNEG,K1,K2,K3,IORD,
     &                       HMAX,ATOL,RTOL,ITOL,RPAR,IPAR)
      NFCN=NFCN+2
      REJECT=.FALSE.
      XOLD=X
      IF (IOUT.GE.1) THEN 
          IRTRN=1 
          HOUT=1.D0
          CALL SOLOUT(NACCPT+1,XOLD,X,Y,N,CONT,ICOMP,NRD,
     &                RPAR,IPAR,IRTRN)
          IF (IRTRN.LT.0) GOTO 79
      END IF
C --- BASIC INTEGRATION STEP  
   1  CONTINUE
      IF (NSTEP.GT.NMAX) GOTO 78
      IF (0.1D0*ABS(H).LE.ABS(X)*UROUND)GOTO 77
      IF ((X+1.01D0*H-XEND)*POSNEG.GT.0.D0) THEN
         H=XEND-X 
         LAST=.TRUE.
      END IF
      NSTEP=NSTEP+1
C --- THE TWELVE STAGES
      IF (IRTRN.GE.2) THEN
         CALL FCN(N,X,Y,K1,RPAR,IPAR)
      END IF
      DO 22 I=1,N 
  22  Y1(I)=Y(I)+H*A21*K1(I)  
      CALL FCN(N,X+C2*H,Y1,K2,RPAR,IPAR)
      DO 23 I=1,N 
  23  Y1(I)=Y(I)+H*(A31*K1(I)+A32*K2(I))  
      CALL FCN(N,X+C3*H,Y1,K3,RPAR,IPAR)
      DO 24 I=1,N 
  24  Y1(I)=Y(I)+H*(A41*K1(I)+A43*K3(I))  
      CALL FCN(N,X+C4*H,Y1,K4,RPAR,IPAR)
      DO 25 I=1,N 
  25  Y1(I)=Y(I)+H*(A51*K1(I)+A53*K3(I)+A54*K4(I))
      CALL FCN(N,X+C5*H,Y1,K5,RPAR,IPAR)
      DO 26 I=1,N 
  26  Y1(I)=Y(I)+H*(A61*K1(I)+A64*K4(I)+A65*K5(I))
      CALL FCN(N,X+C6*H,Y1,K6,RPAR,IPAR)
      DO 27 I=1,N 
  27  Y1(I)=Y(I)+H*(A71*K1(I)+A74*K4(I)+A75*K5(I)+A76*K6(I))
      CALL FCN(N,X+C7*H,Y1,K7,RPAR,IPAR)
      DO 28 I=1,N 
  28  Y1(I)=Y(I)+H*(A81*K1(I)+A84*K4(I)+A85*K5(I)+A86*K6(I)+A87*K7(I))  
      CALL FCN(N,X+C8*H,Y1,K8,RPAR,IPAR)
      DO 29 I=1,N 
  29  Y1(I)=Y(I)+H*(A91*K1(I)+A94*K4(I)+A95*K5(I)+A96*K6(I)+A97*K7(I)
     &   +A98*K8(I))
      CALL FCN(N,X+C9*H,Y1,K9,RPAR,IPAR)
      DO 30 I=1,N 
  30  Y1(I)=Y(I)+H*(A101*K1(I)+A104*K4(I)+A105*K5(I)+A106*K6(I)
     &   +A107*K7(I)+A108*K8(I)+A109*K9(I))
      CALL FCN(N,X+C10*H,Y1,K10,RPAR,IPAR)
      DO 31 I=1,N 
  31  Y1(I)=Y(I)+H*(A111*K1(I)+A114*K4(I)+A115*K5(I)+A116*K6(I)
     &   +A117*K7(I)+A118*K8(I)+A119*K9(I)+A1110*K10(I))
      CALL FCN(N,X+C11*H,Y1,K2,RPAR,IPAR)
      XPH=X+H
      DO 32 I=1,N 
  32  Y1(I)=Y(I)+H*(A121*K1(I)+A124*K4(I)+A125*K5(I)+A126*K6(I)
     &   +A127*K7(I)+A128*K8(I)+A129*K9(I)+A1210*K10(I)+A1211*K2(I))
      CALL FCN(N,XPH,Y1,K3,RPAR,IPAR)
      NFCN=NFCN+11
      DO 35 I=1,N 
      K4(I)=B1*K1(I)+B6*K6(I)+B7*K7(I)+B8*K8(I)+B9*K9(I)
     &   +B10*K10(I)+B11*K2(I)+B12*K3(I)
  35  K5(I)=Y(I)+H*K4(I)
C --- ERROR ESTIMATION  
      ERR=0.D0
      ERR2=0.D0
      IF (ITOL.EQ.0) THEN   
        DO 41 I=1,N 
        SK=ATOLI+RTOLI*MAX(ABS(Y(I)),ABS(K5(I)))
        ERRI=K4(I)-BHH1*K1(I)-BHH2*K9(I)-BHH3*K3(I)
        ERR2=ERR2+(ERRI/SK)**2
        ERRI=ER1*K1(I)+ER6*K6(I)+ER7*K7(I)+ER8*K8(I)+ER9*K9(I)
     &      +ER10*K10(I)+ER11*K2(I)+ER12*K3(I)
  41    ERR=ERR+(ERRI/SK)**2
      ELSE
        DO 42 I=1,N 
        SK=ATOL(I)+RTOL(I)*MAX(ABS(Y(I)),ABS(K5(I)))
        ERRI=K4(I)-BHH1*K1(I)-BHH2*K9(I)-BHH3*K3(I)
        ERR2=ERR2+(ERRI/SK)**2
        ERRI=ER1*K1(I)+ER6*K6(I)+ER7*K7(I)+ER8*K8(I)+ER9*K9(I)
     &      +ER10*K10(I)+ER11*K2(I)+ER12*K3(I)
  42    ERR=ERR+(ERRI/SK)**2
      END IF  
      DENO=ERR+0.01D0*ERR2
      IF (DENO.LE.0.D0) DENO=1.D0
      ERR=ABS(H)*ERR*SQRT(1.D0/(N*DENO))
C --- COMPUTATION OF HNEW
      FAC11=ERR**EXPO1
C --- LUND-STABILIZATION
      FAC=FAC11/FACOLD**BETA
C --- WE REQUIRE  FAC1 <= HNEW/H <= FAC2
      FAC=MAX(FACC2,MIN(FACC1,FAC/SAFE))
      HNEW=H/FAC  
      IF(ERR.LE.1.D0)THEN
C --- STEP IS ACCEPTED  
         FACOLD=MAX(ERR,1.0D-4)
         NACCPT=NACCPT+1
         CALL FCN(N,XPH,K5,K4,RPAR,IPAR)
         NFCN=NFCN+1
C ------- STIFFNESS DETECTION
         IF (MOD(NACCPT,NSTIFF).EQ.0.OR.IASTI.GT.0) THEN
            STNUM=0.D0
            STDEN=0.D0
            DO 64 I=1,N 
               STNUM=STNUM+(K4(I)-K3(I))**2
               STDEN=STDEN+(K5(I)-Y1(I))**2
 64         CONTINUE  
            IF (STDEN.GT.0.D0) HLAMB=ABS(H)*SQRT(STNUM/STDEN) 
            IF (HLAMB.GT.6.1D0) THEN
               NONSTI=0
               IASTI=IASTI+1  
               IF (IASTI.EQ.15) THEN
                  IF (IPRINT.GT.0) WRITE (IPRINT,*) 
     &               ' THE PROBLEM SEEMS TO BECOME STIFF AT X = ',X   
                  IF (IPRINT.LE.0) GOTO 76
               END IF
            ELSE
               NONSTI=NONSTI+1  
               IF (NONSTI.EQ.6) IASTI=0
            END IF
         END IF 
C ------- FINAL PREPARATION FOR DENSE OUTPUT
         IF (IOUT.GE.2) THEN
C ----    SAVE THE FIRST FUNCTION EVALUATIONS   
            DO 62 J=1,NRD
               I=ICOMP(J)
               CONT(J)=Y(I)
               YDIFF=K5(I)-Y(I)
               CONT(J+NRD)=YDIFF
               BSPL=H*K1(I)-YDIFF
               CONT(J+NRD*2)=BSPL
               CONT(J+NRD*3)=YDIFF-H*K4(I)-BSPL
               CONT(J+NRD*4)=D41*K1(I)+D46*K6(I)+D47*K7(I)+D48*K8(I)
     &                  +D49*K9(I)+D410*K10(I)+D411*K2(I)+D412*K3(I)
               CONT(J+NRD*5)=D51*K1(I)+D56*K6(I)+D57*K7(I)+D58*K8(I)
     &                  +D59*K9(I)+D510*K10(I)+D511*K2(I)+D512*K3(I)
               CONT(J+NRD*6)=D61*K1(I)+D66*K6(I)+D67*K7(I)+D68*K8(I)
     &                  +D69*K9(I)+D610*K10(I)+D611*K2(I)+D612*K3(I)
               CONT(J+NRD*7)=D71*K1(I)+D76*K6(I)+D77*K7(I)+D78*K8(I)
     &                  +D79*K9(I)+D710*K10(I)+D711*K2(I)+D712*K3(I)
   62       CONTINUE 
C ---     THE NEXT THREE FUNCTION EVALUATIONS
            DO 51 I=1,N 
  51           Y1(I)=Y(I)+H*(A141*K1(I)+A147*K7(I)+A148*K8(I)
     &            +A149*K9(I)+A1410*K10(I)+A1411*K2(I)+A1412*K3(I)
     &            +A1413*K4(I))
            CALL FCN(N,X+C14*H,Y1,K10,RPAR,IPAR)
            DO 52 I=1,N 
  52           Y1(I)=Y(I)+H*(A151*K1(I)+A156*K6(I)+A157*K7(I)
     &            +A158*K8(I)+A1511*K2(I)+A1512*K3(I)+A1513*K4(I)
     &            +A1514*K10(I))
            CALL FCN(N,X+C15*H,Y1,K2,RPAR,IPAR)
            DO 53 I=1,N 
  53           Y1(I)=Y(I)+H*(A161*K1(I)+A166*K6(I)+A167*K7(I)
     &            +A168*K8(I)+A169*K9(I)+A1613*K4(I)+A1614*K10(I)
     &            +A1615*K2(I))
            CALL FCN(N,X+C16*H,Y1,K3,RPAR,IPAR)
            NFCN=NFCN+3 
C ---     FINAL PREPARATION
            DO 63 J=1,NRD
               I=ICOMP(J)
               CONT(J+NRD*4)=H*(CONT(J+NRD*4)+D413*K4(I)+D414*K10(I)
     &            +D415*K2(I)+D416*K3(I))
               CONT(J+NRD*5)=H*(CONT(J+NRD*5)+D513*K4(I)+D514*K10(I)
     &            +D515*K2(I)+D516*K3(I))
               CONT(J+NRD*6)=H*(CONT(J+NRD*6)+D613*K4(I)+D614*K10(I)
     &            +D615*K2(I)+D616*K3(I))
               CONT(J+NRD*7)=H*(CONT(J+NRD*7)+D713*K4(I)+D714*K10(I)
     &            +D715*K2(I)+D716*K3(I))
  63        CONTINUE
            HOUT=H
         END IF
         DO 67 I=1,N
         K1(I)=K4(I)
  67     Y(I)=K5(I)
         XOLD=X
         X=XPH
         IF (IOUT.GE.1) THEN
            CALL SOLOUT(NACCPT+1,XOLD,X,Y,N,CONT,ICOMP,NRD,
     &                  RPAR,IPAR,IRTRN)
            IF (IRTRN.LT.0) GOTO 79
         END IF 
C ------- NORMAL EXIT
         IF (LAST) THEN
            H=HNEW
            IDID=1
            RETURN
         END IF
         IF(ABS(HNEW).GT.HMAX)HNEW=POSNEG*HMAX  
         IF(REJECT)HNEW=POSNEG*MIN(ABS(HNEW),ABS(H))
         REJECT=.FALSE. 
      ELSE  
C --- STEP IS REJECTED   
         HNEW=H/MIN(FACC1,FAC11/SAFE)
         REJECT=.TRUE.  
         IF(NACCPT.GE.1)NREJCT=NREJCT+1   
         LAST=.FALSE.
      END IF
      H=HNEW
      GOTO 1
C --- FAIL EXIT
  76  CONTINUE
      IDID=-4
      RETURN
  77  CONTINUE
      IF (IPRINT.GT.0) WRITE(IPRINT,979)X   
      IF (IPRINT.GT.0) WRITE(IPRINT,*)' STEP SIZE TOO SMALL, H=',H
      IDID=-3
      RETURN
  78  CONTINUE
      IF (IPRINT.GT.0) WRITE(IPRINT,979)X   
      IF (IPRINT.GT.0) WRITE(IPRINT,*)
     &     ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED' 
      IDID=-2
      RETURN
  79  CONTINUE
      IF (IPRINT.GT.0) WRITE(IPRINT,979)X
 979  FORMAT(' EXIT OF DOP853 AT X=',E18.4) 
      IDID=2
      RETURN
      END
C
      FUNCTION HINIT(N,FCN,X,Y,XEND,POSNEG,F0,F1,Y1,IORD,
     &                       HMAX,ATOL,RTOL,ITOL,RPAR,IPAR)
C ----------------------------------------------------------
C ----  COMPUTATION OF AN INITIAL STEP SIZE GUESS
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),Y1(N),F0(N),F1(N),ATOL(*),RTOL(*)
      DIMENSION RPAR(*),IPAR(*)
C ---- COMPUTE A FIRST GUESS FOR EXPLICIT EULER AS
C ----   H = 0.01 * NORM (Y0) / NORM (F0)
C ---- THE INCREMENT FOR EXPLICIT EULER IS SMALL
C ---- COMPARED TO THE SOLUTION
      DNF=0.0D0
      DNY=0.0D0 
      ATOLI=ATOL(1)
      RTOLI=RTOL(1)    
      IF (ITOL.EQ.0) THEN   
        DO 10 I=1,N 
        SK=ATOLI+RTOLI*ABS(Y(I))
        DNF=DNF+(F0(I)/SK)**2
  10    DNY=DNY+(Y(I)/SK)**2 
      ELSE
        DO 11 I=1,N 
        SK=ATOL(I)+RTOL(I)*ABS(Y(I))
        DNF=DNF+(F0(I)/SK)**2
  11    DNY=DNY+(Y(I)/SK)**2 
      END IF
      IF (DNF.LE.1.D-10.OR.DNY.LE.1.D-10) THEN
         H=1.0D-6
      ELSE
         H=SQRT(DNY/DNF)*0.01D0 
      END IF
      H=MIN(H,HMAX)
      H=SIGN(H,POSNEG) 
C ---- PERFORM AN EXPLICIT EULER STEP
      DO 12 I=1,N
  12  Y1(I)=Y(I)+H*F0(I)
      CALL FCN(N,X+H,Y1,F1,RPAR,IPAR) 
C ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION
      DER2=0.0D0
      IF (ITOL.EQ.0) THEN   
        DO 15 I=1,N 
        SK=ATOLI+RTOLI*ABS(Y(I))
  15    DER2=DER2+((F1(I)-F0(I))/SK)**2   
      ELSE
        DO 16 I=1,N 
        SK=ATOL(I)+RTOL(I)*ABS(Y(I))
  16    DER2=DER2+((F1(I)-F0(I))/SK)**2   
      END IF
      DER2=SQRT(DER2)/H
C ---- STEP SIZE IS COMPUTED SUCH THAT
C ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
      DER12=MAX(ABS(DER2),SQRT(DNF))
      IF (DER12.LE.1.D-15) THEN
         H1=MAX(1.0D-6,ABS(H)*1.0D-3)
      ELSE
         H1=(0.01D0/DER12)**(1.D0/IORD) 
      END IF
      H=MIN(100*H,H1,HMAX)
      HINIT=SIGN(H,POSNEG)  
      RETURN
      END 
C
      FUNCTION CONTD8(II,X,CON,ICOMP,ND)
C ----------------------------------------------------------
C     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT IN CONNECTION
C     WITH THE OUTPUT-SUBROUTINE FOR DOP853. IT PROVIDES AN
C     APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X.
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CON(8*ND),ICOMP(ND)
      COMMON /CONDO8/XOLD,H
C ----- COMPUTE PLACE OF II-TH COMPONENT 
      I=0 
      DO 5 J=1,ND 
      IF (ICOMP(J).EQ.II) I=J
   5  CONTINUE
      IF (I.EQ.0) THEN
         WRITE (6,*) ' NO DENSE OUTPUT AVAILABLE FOR COMP.',II 
         RETURN
      END IF  
      S=(X-XOLD)/H
      S1=1.D0-S
      CONPAR=CON(I+ND*4)+S*(CON(I+ND*5)+S1*(CON(I+ND*6)+S*CON(I+ND*7)))
      CONTD8=CON(I)+S*(CON(I+ND)+S1*(CON(I+ND*2)+S*(CON(I+ND*3)
     &        +S1*CONPAR)))
      RETURN
      END
