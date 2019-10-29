

C     *******************************************************
C     PROGRAM UPMWAKE TO CALCULATE WIND-TURBINE WAKES
C     *******************************************************
C
C       INPUT DATA
C       ==========
C
C     Orography file.............................. OROG$
C     Option 3(0) or 7(1) ec...................... PL
C     Grid size in X direction:................... DX
C     Surface roughness (m):...................... Z0
C     Monin-Obukhov length (m):................... L
C     Ground temperature (K):..................... TS
C     Size of boundary conditions..................NDIAM
C         N�mero de di�metros donde imponer las condiciones
C         de contorno a los lados y arriba del dominio
C     Coeficiente de empuje maximo permitido...... CTMAX
C         Si CT sale mayor se limita a CTMAX y se cambia
C         el diametro del rotor del aerogenerador para
C         que la fuerza no var�a respecto a la que existir�a
C         con el CT y diametros reales
C     Error in Poisson eq. iter.:................. ERRM
C     Error in pressure distrib.:................. ERRORP
C     Relaxation parameter for pressure:.......... RELAJ
C     Maximum number of iterations
C      in Poisson eq. for pressure:............... ITM
C     Number of sections behind the last turbine.. NSL
C     Number of sections to save results?......... NS
C     Output options:
C            SECMAT (NS,1)= NS   Prints results in sec.NS
C            SECMAT (NS,2)= 1    Prints u distrib.
C            SECMAT (NS,3)= 1    Prints v distrib.
C            SECMAT (NS,4)= 1    Prints w distrib.
C            SECMAT (NS,5)= 1    Prints p distrib.
C            SECMAT (NS,6)= 1    Prints t distrib.
C            SECMAT (NS,7)= 1    Prints k distrib.
C            SECMAT (NS,8)= 1    Prints e distrib.
C            SECMAT (NS,9)= 1    Prints turb. visc.
C            SECMAT (NS,10)= 1   Prints module of velocity distrib (horizontal plane) (aun en desarrollo)
C     Number of turbines:..........................S
C  Turbines' data:
C     Coordinate x.................................X1()
C     Coordinate y.................................Y1()
C     Coordinate z.................................Z2()
C     Turbine file.................................TURB$()
C  Directions and incident speeds
C     Meteorological mast height...................HM
C     Number of cases..............................ND
C     Wind direction(clockwise):...................DI()
C     Wind incident speed (In meteor. mast):.......VH()
C  Turbulence Spectra
C     Number of points to be analysed..............NTS
C     Coordinate x (m).............................XTS()
C     Coordinate y (m).............................YTS()
C     Coordinate z (m) ............................ZTS()
C  Wake attenuation
C     Lenght attenuation (number of turbine diameters).......TDWA
C  Meandering calculation
C     Meandering calculation.......................CASO: yes (1) no (0)
C     Meandering type..............................KIND: sinusoidal (1) or spectrum (2)
C     Meandering angle (�) (KIND=1)................TETA0
C     Number of subangles (KIND=1).................NANG
C     Number of time air cross the windpark........PER
C     Time interval definition.....................TID: dt=TID*DIAMETER/VELOCITY
C     Spectrum type (1=Kaimal, 2=Von Karman).......ST
C     Integral Scale Paramenter x direction (m)....LSX
C     Variance x direction (m2/s2).................VA
C     Lower frequency (Hz).........................LF
C     Upper frequency (Hz).........................UF
C     Number of frequency..........................NF
C     Average type.................................ponde: last air cross (0) all cases (1)
C     Calculation of Kinematic Model For Wake Meandering calculation....KMWM: yes (1) no (0)
C            Number of sections to calculate KMWM..NSKMWM
C            Section to calculate KMWM.............SECKMWM



C    ****************************************
C           COMIENZO DEL PROGRAMA
C    ****************************************


C DECLARACION DE VARIABLES

        IMPLICIT NONE

C PARAMETROS

        REAL*8 DENS, BETA, G, PR
        PARAMETER (DENS=1.25, BETA= 3.5E-3, G=9.8, PR=1.)

        REAL*8 SK, SE, CE1, CE2, CU
        PARAMETER (SK=1., SE=1.3, CE1=1.21, CE2=1.92, CU=0.033)

        REAL*8 COMP !, CTMAX
        PARAMETER (COMP=1.E-5 )!, CTMAX=0.85)

        INTEGER*4 SMAX, NDMAX, NTMAX, DIRMAX, NSMAX, NYMAX, NZMAX,
     $            NTSMAX, NXMAX, RANDOMMAX
        PARAMETER (SMAX=150, NDMAX=36, NTMAX=10, DIRMAX=36, NSMAX=150,
     $             NYMAX=800, NZMAX= 120, NTSMAX=20, NXMAX=400,
     $             RANDOMMAX=10000)

        logical*4  od,ex

C DICCIONARIO DE PARAMETROS UTILES PARA EL USUARIO

C    SMAX -----> N� maximo de turbinas en un parque.
C    NDMAX ----> N� maximo de direcciones en las que se puede analizar
C                un parque cualquiera.
C    NTMAX ----> N� maximo de mastiles meteorologicos que se pueden
C                especificar en el fichero de orografia.
C    DIRMAX ---> N� maximo de sectores angulares que se pueden
C                considerar a la hora de construir el fichero de
C                orografia.
C    NSMAX ----> N� maximo de secciones en las que se puede analizar
C                alguna de las variables que proporciona el metodo K-E.
C    NYMAX, NZMAX ---> N� maximo de nudos que pueden formar la malla que
C                      cubre el frente del parque, en su dimension
C                      horizontal y vertical respectivamente.
C    NTSMAX ---> N� maximo de puntos en los que se puede analizar el
C                espectro de la turbulencia.
C    NXMAX ----> N� m�ximo n�mero de planos perpendiculares a la
C                direcci�n de viento que se pueden estudiar
C    CTMAX ----> Coeficiente de empuje maximo permitido. Si CT sale mayor
C                se limita a CTMAX y se cambia el diametro del rotor del
C                aerogenerador para que la fuerza no var�a respecto a la
C                que existir�a con el CT y diametros reales
C    MODIFICAR-> Almacena los nodos donde se interceptan turbina y los
C                l�mites donde no dejar desarrollar la estela de forma
C                libre para cada turbina (IINF,ISUP,HY,HZ,JINF,JSUP,GR)
C    CONTROLKMWM-> En el caso que se quiera hacer el c�lculo del KINEMATIC MODEL FOR WAKE MEANDERING es necesario hacer tres pasadas en el c�lculo:
C                  * 1: sin meandering para hacer el c�lculo del m�ximo decremento de velocidad en toda la secci�n
C                  * 2: con meandering para calcular el valor medio de decremento de velocidad en cada punto
C                  * 3: con meandering para calcular el valor instantaneo de decremento de velocidad en cada punto
C                       y con el hacer el c�lculo del incremento de energ�a cin�tica turbulenta en cada punto rest�ndole el valor medio en cada punto

C VARIABLES LOGICAS

        LOGICAL*4 OPE

C VARIABLES REALES

        REAL*4 randomreal(1:RANDOMMAX)

        REAL*8 DX, DY, DZ, Z0, L, TS, ERRM, ERRORP, RELAJ, HM, HT, AB1,
     $         VH1, ALT, YMINMA, DISTAN, VM1, AK, TK, CT, MINX, CTMAX,
     $         NDIAM, maxpp, OMEGA, DIAMMINIMO, LF, UF,
     $         DIAMAUXILIAR, TIEMPO, TETA, DTIEMPO, Z2MEDIA , VA,
     $         FREQUENCY, DELTAFREQUENCY,TID,tiempolimite,
     $         DIAMMAXIMO, LSX, LSY, VAX, VAY, TETA0,
     $         SPECTRUMX, SPECTRUMY, VHX, VHY, TETA1, TDWA, TIEMPOKMWM

        REAL*8 X1(SMAX), Y1(SMAX), Z2(SMAX), X2(SMAX), Y2(SMAX),
     $         X3(SMAX+NTSMAX), Y3(SMAX),PCU(SMAX),VM(SMAX), VMAX(SMAX),
     $         DIAM(SMAX), DIAM1(SMAX), DIREC(DIRMAX), VO(SMAX, DIRMAX),
     $         MAST(NTMAX, DIRMAX), DI(NDMAX), VH(NDMAX), VMED(DIRMAX),
     $         VH2(NDMAX), VHMEMORY (NDMAX), DELTAV0 (NSMAX)

        REAL*8 X(NYMAX, NZMAX)

        REAL*8 U(NYMAX, NZMAX), V(NYMAX, NZMAX), W(NYMAX, NZMAX),
     $         P(NYMAX, NZMAX), T(NYMAX, NZMAX), E(NYMAX, NZMAX),
     $         K(NYMAX, NZMAX), M(NYMAX, NZMAX), UT(NYMAX,NZMAX,NXMAX),
     $         KT(NYMAX, NZMAX, NXMAX), ET(NYMAX, NZMAX, NXMAX)

        REAL*8 U1(NYMAX, NZMAX), V1(NYMAX, NZMAX), W1(NYMAX, NZMAX),
     $         T1(NYMAX, NZMAX), E1(NYMAX, NZMAX), K1(NYMAX, NZMAX)

        REAL*8 U2(NYMAX, NZMAX), V2(NYMAX, NZMAX), W2(NYMAX, NZMAX),
     $         T2(NYMAX, NZMAX), E2(NYMAX, NZMAX), K2(NYMAX, NZMAX)

        REAL*8 VI(NYMAX, NZMAX), PP(NYMAX, NZMAX),
     $         DELTAVM(NSMAX, NYMAX, NZMAX),
     $         DELTAK(NSMAX, NYMAX, NZMAX)

        REAL*8 U0(NZMAX), T0(NZMAX), E0(NZMAX), VI0(NZMAX), Z(NZMAX),
     $         DU0(NZMAX), DT0(NZMAX), DK0(NZMAX), DE0(NZMAX),
     $         DVI0(NZMAX), D2U0(NZMAX), D2E0(NZMAX), K0(NZMAX),
     $         D2K0(NZMAX), D2T0(NZMAX)

        REAL*8 A(NYMAX - 1), B(NYMAX - 1), C(NYMAX - 1), D(NYMAX - 1)

        REAL*8 XTS(NTSMAX), YTS(NTSMAX), ZTS(NTSMAX), XTS1(NTSMAX),
     $         YTS1(NTSMAX)

        EXTERNAL nSu

C VARIABLES ENTERAS


        INTEGER*4 PL, ITM, NSL, NS, S, ND, NT, DIR1, NM1, NY, NZ, NP,
     $            FLAG1, NS1, HY, HZ, N, I, J, A1, SEC, R, NTS, NANG,
     $            I2, J2, kk, CASO, ctiempo, CASOSTIEMPO, ponde, KIND,
     $            I10, J10, K10, PER, randominteger, NF, INDEXF, ST,
     $            JINF, JSUP, GR, KMWM, NSKMWM, CONTROLKMWM, IJ

        INTEGER*4 SECMAT(NSMAX,10), FLAG(SMAX), FICH(SMAX), DIR(NDMAX),
     $            NM(DIRMAX), NS2(10), WARD(SMAX), SECKMWM(NSMAX)

        INTEGER*4, dimension(:,:), allocatable :: MODIFICAR

C VARIABLES TIPO CARACTER


        CHARACTER*40 OROG$, DATA$, PARQUE$, BSFLOW$,
     $               GRU$, GRV$, GRW$, GRP$, GRK$, GRE$, GRT$, GVI$,
     $               SPC$, IAD$, CTCALCULO$, GRM$, SAL$

        CHARACTER*16 POT$(SMAX),TURB$(SMAX)

        CHARACTER*17 fichero17

        CHARACTER*18 fichero18

        CHARACTER*9 AUXILIAR9

        CHARACTER*80 A$


        CONTROLKMWM=0
        DELTAV0=0.
        DELTAVM=0.
        DELTAK=0.
        secmat=0

        open(178,file='random.dat')
        I=0
        do randominteger=1,RANDOMMAX,1
           I=I+1
           randomreal(I)=360.*RAN(randominteger)
           write(178,*)randomreal(I)
        enddo
        endfile(178)
        close(178)

C ***** INTRODUCCION DEL NOMBRE DEL FICHERO DE DATOS DEL PARQUE *****

1       WRITE (*, *) 'PARK DATA FILE? '

        OPEN(92,FILE='auxiliar.dat')
        READ (92,'(A)', ERR=3) DATA$
        WRITE(*,*) DATA$
        GO TO 2
3       READ(*,'(A)', ERR=1) DATA$

2       R = INDEX (DATA$,'.dat')

        IF (R.EQ.0) THEN
           WRITE (*,*) 'REMEMBER ".dat"'
           GOTO 1
        END IF

        R = R-1

C ****** LECTURA FICHERO INICIALIZACION DEL PARQUE *****


 1111   OPEN (1, FILE= DATA$)

        READ (1, '(A12)') OROG$
        READ (1, '(I3)') PL
        READ (1, '(F8.3)') DX
        READ (1, '(F8.6)') Z0
        READ (1, '(F8.3)') L
        READ (1, '(F8.3)') TS
        READ (1, '(F8.3)') NDIAM
        READ (1, '(F8.6)') CTMAX
        READ (1, '(F8.3)') ERRM
        READ (1, '(F8.3)') ERRORP
        READ (1, '(F8.3)') RELAJ
        READ (1, '(I3)') ITM
        READ (1, '(I3)') NSL
        READ (1, '(I3)') NS
        READ (1, '(A)') A$

        IF ((NS.NE.0).OR.(NS.NE.(-1))) THEN
           DO SEC = 1 , NS
              READ (1, *) (SECMAT(SEC, J), J = 1 , 9)
           END DO
        END IF
        IF (NS.EQ.(-1)) THEN
           READ (1, *) (SECMAT(1, J), J = 1 , 9)
        ENDIF

        READ (1, '(I3)') S
        Z2MEDIA=0.
        DO I = 1 , S
           READ (1, '(F9.1)') X1(I)
           READ (1, '(F9.1)') Y1(I)
           READ (1, '(F8.3)') Z2(I)
           Z2MEDIA=Z2MEDIA+Z2(I)/S
           READ (1, '(A16)') TURB$(I)
        END DO

        READ (1, '(F8.3)') HM
        READ (1, '(I3)') ND
        DO I = 1 , ND
           READ (1, '(F8.3)') DI(I)
           READ (1, '(F8.3)') VH(I)
           VHMEMORY(I)=VH(I)
        END DO

        READ (1,'(I3)') NTS
        DO I = 1 , NTS
           READ (1, '(F8.3)') XTS(I)
           READ (1, '(F8.3)') YTS(I)
           READ (1, '(F8.3)') ZTS(I)
        END DO
        READ (1, '(F8.3)') TDWA
        READ (1, '(I3)') CASO
        OPEN (2 ,FILE=turb$(i))
        READ (2 ,'(F8.3)')DIAMMINIMO
        DIAMMAXIMO=DIAMMINIMO
        CLOSE(2)
        IF (S.GE.2) THEN
           DO I=2,S
              OPEN (2 , FILE=turb$(i))
              READ (2 ,'(F8.3)')DIAMAUXILIAR
              DIAMMINIMO=MIN(DIAMMINIMO,DIAMAUXILIAR)
              DIAMMAXIMO=MAX(DIAMMAXIMO,DIAMAUXILIAR)
              CLOSE(2)
           ENDDO
        ENDIF
        TETA=0.
        IF (CASO.EQ.1) THEN
           IF (PL.EQ.1) THEN
              WRITE(*,fmt='(2a)')'Programa no preparado para ',
     +                     'calcular meandering con 7 ecuaciones'
              STOP
           ENDIF
           READ (1, '(I3)') KIND
           IF ((KIND.NE.1).AND.(KIND.NE.2)) THEN
           WRITE (*,*) 'NO SE SELECCIONA UN TIPO CORRECTO DE MEANDERING'
              STOP
           ENDIF
           TETA0=0.
           READ (1, '(F8.3)') TETA0
!	     TETA0=(180.*ASIN(1.9/(2.5*LOG(Z2MEDIA/Z0))))/3.141592654
           NANG=1
           READ (1,'(I3)') NANG
           READ (1,'(I3)') PER
           READ (1 ,'(F8.3)') TID
           READ (1,'(I3)') ST
           IF ((ST.NE.1).AND.(ST.NE.2)) THEN
              WRITE (*,*) 'NO SE SELECCIONA UN TIPO CORRECTO DE SPECTRO'
              STOP
           ENDIF
           READ (1 ,'(F8.3)') LSX
           IF (ST.EQ.1) THEN
              LSY=2.7*LSX/8.1
              ELSE
                 LSY=LSX
           ENDIF
           READ (1 ,'(F8.3)') VA
           IF (ST.EQ.1) THEN
              VAX=VA
              VAY=0.8*0.8*VA
              ELSE
                 VAX=VA
                 VAY=VA
           ENDIF
           READ (1 ,'(F9.4)') LF
           READ (1 ,'(F9.4)') UF
           DO I=1,ND,1
              IF (UF.GT.(VH(I)/DIAMMAXIMO)) THEN
                 WRITE(*,FMT='(3A,F8.3,A,F8.3,A)')'EL L�MITE SUPERIOR',
     +           ' DE FRECUENCIAS DE C�LCULO CONSIDERADO HASTA ',
     +           'EL MOMENTO, ',UF,' Hz, PARA LA VELOCIDAD ',VH(I),
     +           ', ES MAYOR QUE LA CALCULADA POR EL FILTRADO (fD/U<1)' ! (fLx/U<1)'
                 UF=VH(I)/DIAMMAXIMO !UF=VH(I)/LSX
                 WRITE(*,FMT='(3A,F8.3,A)')'POR TANTO, EL L�MITE ',
     +           'SUPERIOR DE FRECUENCAS DE C�LCULO PARA TODOS LOS',
     +           ' CASOS SE LIMITA A ', UF,' Hz'
!                 pause
              ENDIF
           ENDDO
           READ (1,'(I3)') NF
           READ (1,'(I3)') ponde
           IF ((ponde.NE.0).AND.(ponde.NE.1)) THEN
              WRITE (*,*) 'TIPO DE PONDERACI�N INVALIDO'
              STOP
           ENDIF
           READ (1,'(I3)') KMWM
           IF ((KMWM.NE.0).AND.(KMWM.NE.1)) THEN
              WRITE (*,*) 'ENTRADA DE KMWM NO V�LIDA'
              STOP
           ENDIF
           IF (KMWM.EQ.1) THEN
              CONTROLKMWM=CONTROLKMWM+1
              TIEMPOKMWM=0.
              SELECT CASE (CONTROLKMWM)
                 CASE (1)
                    WRITE(*,FMT='(3a)')" ESTUDIO DEL CASO 1 DEL KMWM:
     +sin meandering para hacer el c�lculo del m�ximo decremento de
     +velocidad en toda la secci�n"
!                    PAUSE
                 CASE (2)
                    WRITE(*,FMT='(3a)')" ESTUDIO DEL CASO 2 DEL KMWM:
     +con meandering para calcular el valor medio de decremento de
     +velocidad en cada punto"
!                    PAUSE
                 CASE (3)
                    WRITE(*,FMT='(5a)')" ESTUDIO DEL CASO 3 DEL KMWM:
     +con meandering para calcular el valor instantaneo de decremento
     +de velocidad en cada punto y con el hacer el c�lculo del
     +incremento de energ�a cin�tica turbulenta en cada punto
     +rest�ndole el valor medio en cada punto"
!                    PAUSE
              END SELECT
              READ (1,'(I3)') NSKMWM
              DO I=1,NSKMWM,1
                 READ (1,'(I3)') SECKMWM(I)
              ENDDO
           ENDIF
           OPEN(93,FILE='INDICE.DAT')
           WRITE (93,'(A)')DATA$
           IF (KIND.EQ.2) THEN
              OPEN(95, FILE='CASOSESPECTRO.DAT')
              WRITE(95,fmt='(2a)')' VINICIAL SPECTRUMX SPECTRUMY ',
     +           'DELTAFREQUENCY TIEMPO VX VY V TETA'
           ENDIF
        ENDIF

        CLOSE (1)

        IF (CONTROLKMWM.EQ.1) CASO=0

C  ****************************
C  *** BUCLE DE DIRECCIONES ***
C  ****************************

        A1 = 0

1004    A1 = A1 + 1

        AB1 = DI(A1)

        CASOSTIEMPO=-1

        IF (KIND.EQ.1) THEN
           OMEGA=0.5*VH(A1)/DIAMMINIMO
!           OMEGA=2.*3.141592654*0.5*VH(A1)/DIAMMINIMO
           DTIEMPO=2.*3.141592654/(OMEGA*NANG)
           ELSEIF (KIND.EQ.2) THEN
              DTIEMPO=TID*DIAMMINIMO/VH(A1)
        ENDIF
        TIEMPO=-1.*DTIEMPO
!        DO J10=1,NYMAX,1
!           DO K10=1,NZMAX,1
!              DO I10=1,NXMAX,1
!                 UT(J10, K10, I10)=0.
!                 KT(J10, K10, I10)=0.
!                 ET(J10, K10, I10)=0.
!              ENDDO
!           ENDDO
!        ENDDO

1005    TIEMPO=TIEMPO+DTIEMPO
        ctiempo=TIEMPO*100.
        CASOSTIEMPO=CASOSTIEMPO+1
        TIEMPOKMWM=TIEMPOKMWM+DTIEMPO

        randominteger=0

C  *********************************************************************************************
C  *** MODIFICACION DE LA VELOCIDAD Y DE LA DIRECCI�N PARA EL CASO DE ESTUDIO DEL MEANDERING ***
C  *********************************************************************************************

        IF (CASO.EQ.1) THEN

           IF (KIND.EQ.1) TETA=TETA0*SIN(OMEGA*TIEMPO)

           IF (KIND.EQ.2) THEN
           VHX=VHMEMORY(A1)
           VHY=0.
           DO INDEXF=1,NF-1,1
              DELTAFREQUENCY=
     +         10.**(LOG10(LF)+(INDEXF)*(LOG10(UF)-LOG10(LF))/(NF-1))-
     +         10.**(LOG10(LF)+(INDEXF-1)*(LOG10(UF)-LOG10(LF))/(NF-1))
              FREQUENCY=0.5*
     +        (10.**(LOG10(LF)+(INDEXF)*(LOG10(UF)-LOG10(LF))/(NF-1))+
     +         10.**(LOG10(LF)+(INDEXF-1)*(LOG10(UF)-LOG10(LF))/(NF-1)))
              IF (ST.EQ.1) THEN
                 SPECTRUMX=4.*VAX*LSX/(VHMEMORY(A1)*
     +              (1+6.*FREQUENCY*LSX/VHMEMORY(A1))**(5./3.))
                 SPECTRUMY=4.*VAY*LSY/(VHMEMORY(A1)*
     +              (1+6.*FREQUENCY*LSY/VHMEMORY(A1))**(5./3.))
                 ELSE
                    SPECTRUMX=4.*VAX*LSX/(VHMEMORY(A1)*
     +              (1+70.8*(FREQUENCY*LSX/VHMEMORY(A1))**2.)**(5./6.))
                    SPECTRUMY=
     +              2.*VAY*LSY*(1+189.*(FREQUENCY*LSY/VHMEMORY(A1))**2.)
     +              /(VHMEMORY(A1)*
     +              (1+70.8*(FREQUENCY*LSY/VHMEMORY(A1))**2.)**(11./6.))
              ENDIF
              randominteger=randominteger+1
              IF (RANDOMINTEGER.GT.RANDOMMAX) randominteger=1
              OMEGA=FREQUENCY*2.*3.141592654
              TETA1=OMEGA*TIEMPO+randomreal(randominteger)

              VHX=VHX+COS(TETA1*3.141592654/180)
     +                      *(SPECTRUMX*DELTAFREQUENCY)**0.5
              VHY=VHY+COS(TETA1*3.141592654/180)
     +                      *(SPECTRUMY*DELTAFREQUENCY)**0.5
           ENDDO
           VH(A1)=(VHX**2.+VHY**2.)**0.5
           TETA=(ATAN(VHY/VHX))*180./3.141592654
           WRITE(95,fmt='(9f12.5)')VHMEMORY(A1),SPECTRUMX,SPECTRUMY,
     +          DELTAFREQUENCY,TIEMPO,VHX,VHY,VH(A1),TETA
           ENDIF

       ENDIF

C ***** APERTURA DE FICHEROS DE SALIDA ******

        IF (CASO.EQ.1) THEN

           WRITE(AUXILIAR9,FMT='(a1,I6.6,a2)')'t',ctiempo,'cs'

           BSFLOW$ = DATA$(:R) // AUXILIAR9 // '.bsf'
           PARQUE$ = DATA$(:R) // AUXILIAR9 // '.pow'
           GRU$ = DATA$(:R) // AUXILIAR9 // '.gru'
           GRV$ = DATA$(:R) // AUXILIAR9 // '.grv'
           GRW$ = DATA$(:R) // AUXILIAR9 // '.grw'
           GRP$ = DATA$(:R) // AUXILIAR9 // '.grp'
           GRT$ = DATA$(:R) // AUXILIAR9 // '.grt'
           GRK$ = DATA$(:R) // AUXILIAR9 // '.grk'
           GRE$ = DATA$(:R) // AUXILIAR9 // '.gre'
           GVI$ = DATA$(:R) // AUXILIAR9 // '.gvi'
           SPC$ = DATA$(:R) // AUXILIAR9 // '.spc'
           IAD$ = DATA$(:R) // AUXILIAR9 // '.iad'
           CTCALCULO$ = DATA$(:R) // AUXILIAR9 // '.ctc'
           GRM$ = DATA$(:R) // AUXILIAR9 // '.grm'
           SAL$ = DATA$(:R) // AUXILIAR9 // '.dat'

           ELSE

              BSFLOW$ = DATA$(:R) // '.bsf'
              PARQUE$ = DATA$(:R) // '.pow'
              GRU$ = DATA$(:R) // '.gru'
              GRV$ = DATA$(:R) // '.grv'
              GRW$ = DATA$(:R) // '.grw'
              GRP$ = DATA$(:R) // '.grp'
              GRT$ = DATA$(:R) // '.grt'
              GRK$ = DATA$(:R) // '.grk'
              GRE$ = DATA$(:R) // '.gre'
              GVI$ = DATA$(:R) // '.gvi'
              SPC$ = DATA$(:R) // '.spc'
              IAD$ = DATA$(:R) // '.iad'
              CTCALCULO$ = DATA$(:R) // '.ctc'
              GRM$ = DATA$(:R) // '.grm'
              SAL$ = DATA$(:R) // '.dat'

        ENDIF

!        BSFLOW$ (14:16)= 'oto'
!        do I=14,24
!           PARQUE$ (I:I) = PARQUE$(I+4:I+4)
!        enddo
!        PARQUE$(25:28)='    '
!        GRU$ (14:16)= 'oto'
!        GRV$(14:16)= 'oto'
!        GRW$(14:16)= 'oto'
!        GRP$(14:16)= 'oto'
!        GRT$(14:16)= 'oto'
!        GRK$(14:16)= 'oto'
!        GRE$(14:16)= 'oto'
!        GVI$(14:16)= 'oto'
!        SPC$(14:16)= 'oto'
!        IAD$(14:16)= 'oto'
!        CTCALCULO$(14:16)= 'oto'
!        GRM$ (14:16)= 'oto'
!        SAL$ (14:16)= 'oto'

C ***** DICCIONARIO DE NUMEROS ASOCIADOS A CADA FICHERO *****

C       DATA$ -------> 1
C       BSFLOW$ -----> 2
C       PARQUE$ -----> 3
C       PANTALLA ----> 6
C       GRU$ --------> 7
C       GRV$ --------> 8
C       GRW$ --------> 9
C       GRP$ --------> 10
C       GRT$ --------> 11
C       GRK$ --------> 12
C       GRE$ --------> 13
C       GVI$ --------> 14
C       SPC$ --------> 15
C       OROG$ -------> 16
C       GRM$ --------> 17
C       TURB$ -------> 19 en adelante, tantos como turbinas.
C       IAD$ --------> 90
C       CTCALCULO$ --> 91
C       FICHERO DE DATOS -> 92
C       Indice de casos ->93
C       SAL$ --------> 94
C       Casos espectro -> 95

!        OPEN (94, FILE= SAL$)

!        WRITE (94, '(A12)') OROG$
!        WRITE (94, '(I3)') PL
!        WRITE (94, '(F8.3)') DX
!        WRITE (94, '(F8.6)') Z0
!        WRITE (94, '(G8.3)') L
!        WRITE (94, '(F8.3)') TS
!        WRITE (94, '(F8.3)') NDIAM
!        WRITE (94, '(F8.6)') CTMAX
!        WRITE (94, '(F8.3)') ERRM
!        WRITE (94, '(F8.3)') ERRORP
!        WRITE (94, '(F8.3)') RELAJ
!        WRITE (94, '(I3)') ITM
!        WRITE (94, '(I3)') NSL
!        WRITE (94, '(I3)') NS
!        WRITE (94, '(A)') A$

!        IF ((NS.NE.0).OR.(NS.NE.(-1))) THEN
!           DO SEC = 1 , NS
!              WRITE (94, *) (SECMAT(SEC, J), J = 1 , 9)
!           END DO
!        END IF
!        IF (NS.EQ.(-1)) THEN
!           WRITE (94, *) (SECMAT(1, J), J = 1 , 9)
!        ENDIF

!        WRITE (94, '(I3)') S

!        DO I = 1 , S
!           WRITE (94, '(F9.1)') X1(I)
!           WRITE (94, '(F9.1)') Y1(I)
!           WRITE (94, '(F8.3)') Z2(I)
!           WRITE (94, '(A16)') TURB$(I)
!        END DO

!        WRITE (94, '(F8.3)') HM
!        WRITE (94, '(I3)') ND
!        DO I = 1 , ND
!           WRITE (94, '(F8.3)') DI(I)
!           WRITE (94, '(F8.3)') VH(I)
!        END DO

!        WRITE (94,'(I3)') NTS
!        DO I = 1 , NTS
!           WRITE (94, '(F8.3)') XTS(I)
!           WRITE (94, '(F8.3)') YTS(I)
!           WRITE (94, '(F8.3)') ZTS(I)
!        END DO
!        WRITE (94,'(F8.3)') TDWA
!        WRITE (94, '(I3)') CASO
!        IF (CASO.EQ.1) THEN
!           WRITE (94, '(I3)') KIND
!           WRITE (94, '(F8.3)') TETA0
!           WRITE (94, '(I3)') NANG
!           WRITE (94, '(I3)') PER
!           WRITE (94, '(F8.3)') TID
!           WRITE (94, '(I3)') ST
!           WRITE (94, '(F8.3)') LSX
!           WRITE (94, '(F8.3)') VA
!           WRITE (94,'(F9.4)') LF
!           WRITE (94,'(F9.4)') UF
!           WRITE (94,'(I3)') NF
!           WRITE (94,'(I3)') ponde
!           WRITE (94,'(I3)') KMWM
!           IF (KMWM.EQ.1) THEN
!              WRITE(94,'(I3)') NSKMWM
!              DO I=1,NSKMWM
!                 WRITE(94,'(I3)') SECKMWM(I)
!              ENDDO
!           ENDIF
!        ENDIF

!        ENDFILE(94)
!        CLOSE (94)

        OPEN (2, FILE=BSFLOW$)
        OPEN (3, FILE=PARQUE$)

        DO I = 1 , S
           inquire(file=turb$(i),exist=ex,opened=od)
           IF (.NOT.ex) THEN
              WRITE(*,fmt='(4a)')' Error. El fichero del tipo de ',
     +        'turbina : ',turb$(i),' no existe o no puede ser abierto.'
              stop
           ENDIF
           IF (.NOT.od) then
              FICH(I)=18+I
              OPEN (FICH(I), FILE= TURB$(I))
              ELSE
                 kk=1
		       DO WHILE (TURB$(I).NE.TURB$(kk))
	              kk=kk+1
		       ENDDO
                 FICH(I)=18+kk
           ENDIF
        END DO

        DO J = 2 , 10
          NS2(J) = 0
          DO I = 1 , NS
            IF (SECMAT(I, J).EQ.1) THEN
              NS2(J) = 1
            END IF
          END DO
          IF ((NS.EQ.(-1)).AND.(SECMAT(1, J).EQ.1)) NS2(J) = 1
        END DO

        IF (NS2(2).EQ.1) THEN
            OPEN (7, FILE=GRU$)
        END IF

        IF (NS2(3).EQ.1) THEN
            OPEN (8, FILE=GRV$)
        END IF

        IF (NS2(4).EQ.1) THEN
            OPEN (9, FILE=GRW$)
        END IF

        IF (NS2(5).EQ.1) THEN
            OPEN (10, FILE=GRP$)
        END IF

        IF (NS2(6).EQ.1) THEN
            OPEN (11, FILE=GRT$)
        END IF

        IF (NS2(7).EQ.1) THEN
            OPEN (12, FILE=GRK$)
        END IF

        IF (NS2(8).EQ.1) THEN
            OPEN (13, FILE=GRE$)
        END IF

        IF (NS2(9).EQ.1) THEN
            OPEN (14, FILE=GVI$)
        END IF

        IF (NS2(10).EQ.1) THEN
            OPEN (17, FILE=GRM$)
        END IF

        IF (NTS.NE.0) THEN
            OPEN (15, FILE=SPC$)
        END IF

        OPEN (90, FILE=IAD$)

        OPEN (91, FILE=CTCALCULO$)

C ***** OROGRAFIA: LECTURA DEL FICHERO DE OROGRAFIA *****


        IF (OROG$.NE.'no') THEN

           CALL FICHOROG (OROG$, S, DIR1, NT, VO, DIREC, MAST, VMED,
     $                    DIRMAX, NTMAX)

        END IF


C ***** LECTURA DE DIAMETROS *****

        DO I = 1 , S
          inquire(unit=fich(i),opened=od)
          if(od) then
             REWIND FICH(I)
             READ (FICH(I), '(F8.3)') DIAM(I)
             if (diam(i).eq.0.) then
                write(*,fmt='(a,i2.2,a)')' El di�metro de la turbina ',
     +                  i,' vale cero'
                stop
             endif
             else
                write(*,fmt='(a,i2.2,a)')' Error1: La unidad ',
     +                fich(i),' no est� abierta'
                stop
          end if
        END DO

C ****** CALCULO DE LA VELOCIDAD INCIDENTE CORREGIDA ********


        CALL CORRECVEL (ND, NT, A1, DI, NM, AB1, NM1, DIR1, VH2, VH,
     $                  DIR, OROG$, VH1, DIREC, MAST, NTMAX, DIRMAX)



C ****** ALTURA DEL MASTIL METEOROLOGICO (SI EXISTE) *****


        CALL ALTURA (S, Z2, ALT, HM, HT, COMP)


C ****** ENCABEZAMIENTO PARA CADA CASO DE ANALISIS ******


        CALL ENCABEC (ND, OROG$, DATA$, NM1, ALT, VH, AB1, 2, A1,TIEMPO)
        CALL ENCABEC (ND, OROG$, DATA$, NM1, ALT, VH, AB1, 3, A1,TIEMPO)
        CALL ENCABEC (ND, OROG$, DATA$, NM1, ALT, VH, AB1, 6, A1,TIEMPO)
        CALL ENCABEC (ND, OROG$, DATA$, NM1, ALT, VH, AB1,90, A1,TIEMPO)
        CALL ENCABEC (ND, OROG$, DATA$, NM1, ALT, VH, AB1,91, A1,TIEMPO)

        IF (NS2(2).EQ.1) THEN
        CALL ENCABEC (ND, OROG$, DATA$, NM1, ALT, VH, AB1, 7, A1,TIEMPO)
        END IF

        IF (NS2(3).EQ.1) THEN
        CALL ENCABEC (ND, OROG$, DATA$, NM1, ALT, VH, AB1, 8, A1,TIEMPO)
        END IF

        IF (NS2(4).EQ.1) THEN
        CALL ENCABEC (ND, OROG$, DATA$, NM1, ALT, VH, AB1, 9, A1,TIEMPO)
        END IF

        IF (NS2(5).EQ.1) THEN
        CALL ENCABEC(ND, OROG$, DATA$, NM1, ALT, VH, AB1, 10, A1,TIEMPO)
        END IF

        IF (NS2(6).EQ.1) THEN
        CALL ENCABEC(ND, OROG$, DATA$, NM1, ALT, VH, AB1, 11, A1,TIEMPO)
        END IF

        IF (NS2(7).EQ.1) THEN
        CALL ENCABEC(ND, OROG$, DATA$, NM1, ALT, VH, AB1, 12, A1,TIEMPO)
        END IF

        IF (NS2(8).EQ.1) THEN
        CALL ENCABEC(ND, OROG$, DATA$, NM1, ALT, VH, AB1, 13, A1,TIEMPO)
        END IF

        IF (NS2(9).EQ.1) THEN
        CALL ENCABEC(ND, OROG$, DATA$, NM1, ALT, VH, AB1, 14, A1,TIEMPO)
        END IF

        IF (NS2(10).EQ.1) THEN
        CALL ENCABEC(ND, OROG$, DATA$, NM1, ALT, VH, AB1, 17, A1,TIEMPO)
        END IF

        IF (NTS.NE.0) THEN
        CALL ENCABEC(ND, OROG$, DATA$, NM1, ALT, VH, AB1, 15, A1,TIEMPO)
        END IF


C ****** CAMBIO DE EJES ******


        CALL CAMBIOEJES (S, NTS, X2, Y2, X1, Y1, XTS, YTS, XTS1,
     $                   YTS1, AB1)


C *********************************************
C *******  CREACION DE LA MALLA INICIAL *******
C *********************************************


C ***** CALCULO DE INCREMENTOS DE LA MALLA DY, DZ *****


        CALL INCREMEN (S, DIAM1, DIAM, DY, DZ)

C ***** DETERMINACION DE LOS LIMITES DE LA MALLA *****


        CALL LIMMALLA (X2, Y2, X3, Y3, XTS1, S, NTS, DIAM1, DY, NY, NSL,
     $                 DX, NP, HT, Z0, DZ, NZ, YMINMA, MINX, NDIAM)
        IF (CASO.EQ.1) THEN
           IF (CASOSTIEMPO.EQ.0) THEN
              tiempolimite=NP*DX/VHMEMORY(A1)
              WRITE(93,'(A)')'Nombre de ficheros de casos calculados'
           ENDIF
           WRITE(93,FMT='(A)')DATA$(:R) // AUXILIAR9
        ENDIF


        IF (NS.EQ.(-1)) THEN
           NS = NP
           DO SEC = 1 , NS
              SECMAT(SEC, 1) = SEC
              DO J = 2 , 10
                 SECMAT(SEC, J) = NS2(J)
              END DO
           END DO
        END IF

        IF ((NY.GT.NYMAX).OR.(NZ.GT.NZMAX)) THEN

           WRITE (*, 10001) NY, NZ
10001      FORMAT (/, 1X, 'GRID TOO LARGE',/,1X, 'NY = ', I8, /,  1X,
     $            'NZ = ', I8 )
           STOP

        END IF


C     ****************************
C     ******* FLUJO BASICO *******
C     ****************************



C     CALCULO DEL FLUJO BASICO

        CALL BASICFLOW (NZ, DZ, Z, L, ALT, VH1, Z0, VI0, U0, E0, K0,
     $                  T0, DVI0, DU0, DE0, DK0, DT0, D2U0, D2E0, D2K0,
     $                  D2T0, CU, TS, G, PR, TETA)

C     IMPRESION FICHERO DEL FLUJO BASICO


        CALL IMPBSFLOW (NZ, Z, U0, K0, E0, VI0)


C     *************************************
C     ***** CALCULOS EN CADA SECCION ******
C     *************************************


C     INICIALIZACION DE LOS NUDOS DE LA MALLA


        CALL INICMALL (NY, NZ, U, V, W, P, PP, T, K, E)


C     COMIENZO DEL BUCLE DE SECCIONES


        N = 0
        ALLOCATE(MODIFICAR(1:S,1:7))
        MODIFICAR=-2

!        WRITE (*, 10002) NP
!10002   FORMAT (1X, 'TOTAL NUMBER OF SECTIONS = ', I3)

!   42      WRITE (*, 10003) N
!10003   FORMAT (1X, 'SECTION NUMBER: ', I3)

        DISTAN = N * DX

        CALL CALCULO (PL, NY, NZ, N, ITM, DX, DY, DZ, DENS, BETA, G,
     $         ERRM, ERRORP, RELAJ, PR, SK, SE, CE2, CE1, U, V, W, P,
     $         T, E, K, VI, U0, E0, VI0, DU0, DT0, DK0, DE0, DVI0,
     $         D2U0, D2E0, D2K0, D2T0, K0, U1, V1, W1, T1, E1, K1, U2,
     $         V2, W2, T2, E2, K2, X, A, B, C, D, PP, A1, NP, ND,
     $         DATA$, maxpp, TETA, UT, CASO, DTIEMPO, KT, ET, TIEMPO,
     $         S, MODIFICAR)


        IF (CONTROLKMWM.NE.0) CALL CALCULOKMWM (CONTROLKMWM,NSKMWM,
     +            SECKMWM,N,NY,NZ,DELTAV0,U,DELTAVM,DTIEMPO,DELTAK)

C ***********************************************************
C ********* ANALISIS DEL ESPECTRO DE LA TURBULENCIA *********
C ***********************************************************

        IF (NTS.NE.0) THEN

           CALL ESPECTRO (XTS, YTS, ZTS, XTS1, YTS1, DX, DY, DZ, NTS,
     $                   DISTAN, U, K, E, U0, K0, E0, NY, NZ, YMINMA,
     $                   S, MINX, nSu)

        END IF

C *********************************************
C ********* ANALISIS DE INTERACCION ***********
C *********************************************


        CALL INTERAC (S, MINX, DX, DISTAN, X2, FLAG, FLAG1)

        IF (FLAG1.EQ.0) THEN     ! NO EXISTE INTERACCION
           GO TO 610
        END IF


C **************************************
C ********REDEFINICION DE MALLA*********
C **************************************


        DO I = 1 , S

          IF (FLAG(I).EQ.0) GO TO 600

          WRITE (*, 10004) I
10004     FORMAT (/, 1X, 'INTERACTION WITH TURBINE ', I3)

C       DETERMINACION DE LA POSICION DEL ROTOR DE INTERACCION EN LA MALLA

           HZ = INT(Z2(I) / DZ + .5) + 1
           HY = INT((Y2(I) - YMINMA) / DY + .5) + 1

C       DETERMINACION DE LA VELOCIDAD EN EL CENTRO DEL ROTOR

           CALL VCEN (VM1, U0, U, HY, HZ, NY, NZ)

C       DETERMINACION DE LA K EN EL CENTRO DEL ROTOR

           I2=0
           DO J2 = 1 , S , 1
              IF ((X2(I).EQ.(X1(J2)*COS(AB1*3.141592654/180)-Y1(J2)*
     $             SIN(AB1*3.141592654/180))).AND.(Y2(I).EQ.(Y1(J2)*
     $             COS(AB1*3.141592654/180)+X1(J2)*
     $             SIN(AB1*3.141592654/180)))) I2 = J2
           ENDDO

           CALL KCEN (K, E, HY, HZ, NY, NZ, X1, Y1, I2, S,
     $                Z, U0, K0, E0, HM, Z2)

C       LECTURA DE LAS CARACTERISTICAS DE LA TURBINA INTERCEPTADA

           CALL TURBINA (S, I, FICH, AK, TK, VM1, CT, POT$, WARD,
     $                  CTMAX, DIAM(I))

C       CAMBIO DE LAS CONDICIONES DE CONTORNO EN LA MALLA

           CALL CONTORNO (S, I, NY, NZ, DIAM, CT, DY, HY, HZ, VM, VMAX,
     $                    U0, U, VM1, T, TK, K, AK, E, K0, E0, JINF,
     $                    JSUP, GR)

           MODIFICAR(I,1)=N
           MODIFICAR(I,2)=N+INT(DIAM(I)*TDWA/DX)
           MODIFICAR(I,3)=HY
           MODIFICAR(I,4)=HZ
           MODIFICAR(I,5)=JINF
           MODIFICAR(I,6)=JSUP
           MODIFICAR(I,7)=GR

C       EFECTO DE LA OROGRAFIA EN LA POTENCIA

           IF (OROG$.NE.'no') THEN

              CALL OROGPOT (I, A1, DIR, DIR1, AB1, VM, VMAX, VH1, S,
     $                      VO, DIREC, ND, DIRMAX)

           END IF

600     END DO


C ************************************
C ********CALCULO DE POTENCIA*********
C ************************************


        CALL POTENCIA (S, VM, FLAG, PCU, POT$, FICH, WARD)


C **************************************
C **** ALMACENAMIENTO DE RESULTADOS ****
C **************************************


610     IF (NS.EQ.0) THEN
           GO TO 620
        END IF


C     AVERIGUAMOS SI LA SECCION ACTUAL ESTA ESPECIFICADA EN EL FICHERO
C     DE DEFINICION DEL PARQUE PARA ALMACENAR DATOS

        FLAG1 = 0
C         flag1=1
        DO I = 1 , NS
          IF (SECMAT(I, 1).EQ.N) THEN
            NS1 = I
            FLAG1 = 1
          END IF
        END DO

        IF (FLAG1.EQ.0) THEN
           GO TO 620
        END IF

        IF (SECMAT(NS1, 2).EQ.1) THEN
            CALL RESULTADOS (N, NY, NZ, U, 7)
        END IF

        IF (SECMAT(NS1, 3).EQ.1) THEN
            CALL RESULTADOS (N, NY, NZ, V, 8)
        END IF

        IF (SECMAT(NS1, 4).EQ.1) THEN
            CALL RESULTADOS (N, NY, NZ, W, 9)
        END IF

        IF (SECMAT(NS1, 5).EQ.1) THEN
            CALL RESULTADOS (N, NY, NZ, P, 10)
        END IF

        IF (SECMAT(NS1, 6).EQ.1) THEN
            CALL RESULTADOS (N, NY, NZ, T, 11)
        END IF

        IF (SECMAT(NS1, 7).EQ.1) THEN
            CALL RESULTADOS (N, NY, NZ, K, 12)
        END IF

        IF (SECMAT(NS1, 8).EQ.1) THEN
            CALL RESULTADOS (N, NY, NZ, E, 13)
        END IF

        IF (SECMAT(NS1, 9).EQ.1) THEN
            CALL RESULTADOS (N, NY, NZ, VI, 14)
        END IF

        DO J = 1 , NY
        DO I = 1 , NZ
           M(J, I) = SQRT(U(J,I)**2.+V(J,I)**2.) !�Esta bien? Creo que no ya que esto ser�a el m�dulo de la composici�n de las "perturbaciones" de la velocidad
! No esta implementado, a�n. Habr�a que cambiar las cuatro l�neas donde aparece secmat variando de 1 a 9, que variase de 1 a 10
        END DO
        END DO
        IF (SECMAT(NS1, 10).EQ.1) THEN
            CALL RESULTADOS (N, NY, NZ, M, 17)
        END IF


C ***************************
C *****�ULTIMA SECCION?******
C ***************************


620     IF (N.LT.NP) THEN
           N = N + 1
           GO TO 42
        END IF

        DEALLOCATE(MODIFICAR)

C **********************
C **** RENDIMIENTOS ****
C **********************

        CALL RENDIMIENTOS (PCU, S, VMAX, TURB$, COMP,fich)

        IF (CASO.EQ.1) THEN
           ENDFILE(2)
           CLOSE(2)
           ENDFILE(3)
           CLOSE(3)
           IF (NS2(2).EQ.1) THEN
              ENDFILE(7)
              CLOSE(7)
           ENDIF
           IF (NS2(3).EQ.1) THEN
              ENDFILE(8)
              CLOSE(8)
           ENDIF
           IF (NS2(4).EQ.1) THEN
              ENDFILE(9)
              CLOSE(9)
           ENDIF
           IF (NS2(5).EQ.1) THEN
              ENDFILE(10)
              CLOSE(10)
           ENDIF
           IF (NS2(6).EQ.1) THEN
              ENDFILE(11)
              CLOSE(11)
           ENDIF
           IF (NS2(7).EQ.1) THEN
              ENDFILE(12)
              CLOSE(12)
           ENDIF
           IF (NS2(8).EQ.1) THEN
              ENDFILE(13)
              CLOSE(13)
           ENDIF
           IF (NS2(9).EQ.1) THEN
              ENDFILE(14)
              CLOSE(14)
           ENDIF
           IF (NS2(10).EQ.1) THEN
              ENDFILE(17)
              CLOSE(17)
           ENDIF
           IF (NTS.NE.0) THEN
              ENDFILE(15)
              CLOSE(15)
           ENDIF
           ENDFILE(90)
           CLOSE(90)
           ENDFILE(91)
           CLOSE(91)
        ENDIF

C ********************************
C ****** �ULTIMO INSTANTE?  ******
C ********************************


        IF ((TIEMPO.LT.
     +     (PER*NP*DX/VHMEMORY(A1))).AND.(CASO.EQ.1)) THEN
           GO TO 1005
	  END IF
        IF (CASO.EQ.1) THEN
           ENDFILE(93)
           CLOSE(93)
           IF (KIND.EQ.2) THEN
              ENDFILE(95)
              CLOSE(95)
           ENDIF
        ENDIF

        IF (CONTROLKMWM.EQ.2) DELTAVM=DELTAVM/TIEMPOKMWM
        IF (CONTROLKMWM.EQ.3) DELTAK=DELTAK/TIEMPOKMWM

C ********************************
C ****** �ULTIMA DIRECCION? ******
C ********************************

        IF (A1.LT.ND) THEN
           GO TO 1004
        END IF


c      ===========================================================
c                        FICHERO GRAFICO
c      ===========================================================

        do n=200,208
           if (n.eq.200) then
!        open (unit=n,file='Ficheros\pow\oto\u.header')
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',gru$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csu.header'
!        open (unit=n,file=fichero17)
              open (unit=n,file='u.header')
              write(n,*) 'file = ',gru$
           end if
           if (n.eq.201) then
!        open (unit=n,file='Ficheros\pow\oto\v.header')
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grv$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csv.header'
!        open (unit=n,file=fichero17)
              open (unit=n,file='v.header')
              write(n,*) 'file = ',grv$
           end if
           if (n.eq.202) then
!        open (unit=n,file='Ficheros\pow\oto\w.header')
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grw$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csw.header'
!        open (unit=n,file=fichero17)
              open (unit=n,file='w.header')
              write(n,*) 'file = ',grw$
           end if
           if (n.eq.203) then
!        open (unit=n,file='Ficheros\pow\oto\p.header')
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grp$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csp.header'
!        open (unit=n,file=fichero17)
              open (unit=n,file='p.header')
              write(n,*) 'file = ',grp$
           end if
           if (n.eq.204) then
!        open (unit=n,file='Ficheros\pow\oto\t.header')
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grt$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'cst.header'
!        open (unit=n,file=fichero17)
              open (unit=n,file='t.header')
              write(n,*) 'file = ',grt$
           end if
           if (n.eq.205) then
!        open (unit=n,file='Ficheros\pow\oto\k.header')
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grk$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csk.header'
!        open (unit=n,file=fichero17)
              open (unit=n,file='k.header')
              write(n,*) 'file = ',grk$
           end if
           if (n.eq.206) then
!        open (unit=n,file='Ficheros\pow\oto\e.header')
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',gre$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'cse.header'
!        open (unit=n,file=fichero17)
              open (unit=n,file='e.header')
              write(n,*) 'file = ',gre$
           end if
           if (n.eq.207) then
!        open (unit=n,file='Ficheros\pow\oto\vi.header')
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',gvi$
!        write(fichero18,fmt='(a1,i6.6,a11)')'t',ctiempo,'csvi.header'
!        open (unit=n,file=fichero18)
              open (unit=n,file='vi.header')
              write(n,*) 'file = ',gvi$
           end if
           if (n.eq.208) then
!        open (unit=n,file='Ficheros\pow\oto\m.header')
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grm$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csm.header'
!        open (unit=n,file=fichero17)
              open (unit=n,file='m.header')
              write(n,*) 'file = ',grm$
           end if
           write(n,*) 'grid = ',np,' x ',nz,' x ',ny-2
           write(n,*) 'field = scalar_one'
           write(n,*) 'type = float'
           write(n,*) 'format = ascii'
           write(n,*) 'structure = scalar'
           write(n,*) 'interleaving = record'
           write(n,*) 'header = lines 12'
           write(n,2222)  dx,dz,dy
2222       format('positions = regular, regular, regular,'
     +        ,'0.,',f6.2,',0.,',f6.2,',0.,',f6.2)
           endfile(n)
           close(n)
        end do

C ****************************************************
C ****** PONDERACI�N DE LOS CASOS DE MEANDERING ******
C ****************************************************

        if (caso.eq.1) THEN

           CALL MEANDERING(NP,NY-2,NZ,ponde,per,R,
     $                        dtiempo,tiempolimite,DATA$,S)

        endif


C ********************************
C ******* FIN DEL PROGRAMA *******
C ********************************

        DO I=1,300
           if (i.eq.6) cycle
           INQUIRE (I, OPENED=OPE)
           IF (OPE) CLOSE(I)
        ENDDO

        IF ((1.LE.CONTROLKMWM).AND.(CONTROLKMWM.LE.2)) GOTO 1111

        IF (CONTROLKMWM.NE.0) CALL GUARDARKMWM(NY,NZ,NSKMWM,SECKMWM,
     +                         DY,DZ,U,DELTAV0,DELTAVM,DELTAK)

        STOP
        END


C **************************
C ******* SUBRUTINAS *******
C **************************



        SUBROUTINE ALTURA (S, Z2, ALT, HM, HT, COMP)

C --------------------------------------------------------------------------
C     Subrutina que determina la altura a la que se mide la velocidad in-
C  cidente en el parque en funcion de que exista o no mastil meteorologico.
C --------------------------------------------------------------------------


C       DECLARACION DE VARIABLES

        IMPLICIT NONE

C       VARIABLES INTERNAS

        INTEGER*4 I

C       VARIABLES EXTERNAS

        INTEGER*4 S

        REAL*8 ALT, HM, HT, COMP, Z2(S)

C  *** CALCULO DE LA ALTURA MEDIA DE LAS TURBINAS DEL PARQUE ***

        HT = 0.
        DO I = 1 , S
           HT = HT + Z2(I)
        END DO
        HT = HT / S


C  *** ALTURA DEL MASTIL METEOROLOGICO (SI EXISTE) ***

        IF (ABS(HM).LT.COMP) THEN
           ALT = HT        !SE TOMA COMO ALTURA DEL MASTIL LA ALTURA MEDIA
        ELSE               !DEL ROTOR SI EL PARAMETRO HM (FICHERO DE
           ALT = HM        !DEFINICION DEL PARQUE) ES NULO.
        END IF

        RETURN
        END

C ==========================================================================

        SUBROUTINE BASICFLOW (NZ, DZ, Z, L, ALT, VH1, Z0, VI0, U0, E0,
     $                        K0, T0, DVI0, DU0, DE0, DK0, DT0, D2U0,
     $                        D2E0, D2K0, D2T0, CU, TS, G, PR, TETA)

C  ----------------------------------------------------------------------
C     Subrutina que calcula las magnitudes que describen el flujo basico
C  del lugar en cuestion, para cada caso de analisis.
C  ----------------------------------------------------------------------

C       DECLARACION DE VARIABLES

        IMPLICIT NONE

C       VARIABLES INTERNAS

        INTEGER*4 I, PR0

        REAL*8 X1, FVI, FU, UFR, FE, FK, FT, DFVI, DFE, D2FE, D2FVI

C       VARIABLES EXTERNAS

        INTEGER*4 NZ

        REAL*8 DZ, L, ALT, VH1, Z0, CU, TS, G, PR, TETA

        REAL*8 U0(NZ), T0(NZ), E0(NZ), VI0(NZ), Z(NZ), DU0(NZ), DT0(NZ),
     $         DK0(NZ), DE0(NZ), DVI0(NZ), D2U0(NZ), D2E0(NZ), K0(NZ),
     $         D2K0(NZ), D2T0(NZ)

C   COMIENZO

        Z(1) = Z0
        DO I = 2 , NZ
           Z(I) = DZ * (I - 1)
        END DO

C   IMPOSICION DE CONDICION DE LA VELOCIDAD INCIDENTE. CALCULO DE LA
C   VELOCIDAD DE FRICCION.

        IF (L.LE.0.) THEN
           FVI = (1. - 16. * ALT / L) ** (-.25)
           X1 = 1. / FVI
           FU = LOG((1. + X1 ** 2) / 2. * ((1. + X1) / 2.) ** 2)-2.*
     $          ATAN(X1) + ATAN(1.) * 2.
        ELSE
           FU = -5. * ALT / L
        END IF

               !IMPOSICION DE LA CONDICION
               !DE CONTORNO V(ALT)=VH1 --> UFR

        UFR = VH1 / (2.5 * (LOG(ALT / Z0) - FU))

        DO I = 1 , NZ

        IF (L.LE.0.) THEN
          PR0=2
          FVI = (1 - 16 * Z(I) / L) ** (-.25)
          X1 = 1. / FVI
          FU = LOG((1. + X1 ** 2) / 2. * ((1. + X1) / 2.) ** 2)
     $                - 2. * ATAN(X1) + ATAN(1.) * 2
          FE = 1. - Z(I) / L
          FK = (FE / FVI) ** .5
          FT = 2. * LOG((1. + X1 ** 2) / 2.)
          DFVI = .25 * (1. - 16. * Z(I) / L) ** (-1.25) * 16. / L
          D2FVI = .3125 * (16./L)**2 * (1-16.*Z(I)/L)**(-2.25)
          DFE = -1. / L
          D2FE = 0.
        ELSE
          PR0=1
          FVI = 1. + 5. * Z(I) / L
          FU = -5. * Z(I) / L
          FE = (1. + 2.5 * (Z(I) / L) ** .6) ** 1.5
          FK = (FE / FVI) ** .5
          FT = -5. * Z(I) / L
          DFVI = 5. / L
          D2FVI = 0.
          DFE = FE ** (1. / 3.) * 2.5 / L * (Z(I) / L) ** (-.4)
          D2FE = 2.5 / L * (-.4 / L ** (-.4) * (Z(I)) ** (-1.4)
     $           * FE ** (1. / 3.) + (Z(I) / L) ** (-.4) * FE ** (2./3.)
     $           * DFE / 3)
        END IF

C   FLUJO BASICO: CALCULO DE VARIABLES Y SUS DERIVADAS

        VI0(I) = .4 * UFR * Z(I) / FVI
        U0(I) = 2.5 * UFR * (LOG(Z(I) / Z0) - FU) *
     $                COS (TETA * 3.141592654 / 180)
        E0(I) = UFR ** 3 / .4 / Z(I) * FE
        K0(I) = UFR ** 2 * FK / CU ** .5
        T0(I) = TS + 2.5 * (LOG(Z(I) / Z0) - FT) * UFR ** 2
     $          * TS / L / .4 / G
        DVI0(I) = .4 * UFR * (FVI - Z(I) * DFVI) / FVI ** 2
        DU0(I) = 2.5 * UFR / Z(I) * FVI
        DE0(I) = 2.5 * UFR ** 3 * (DFE * Z(I) - FE) / Z(I) ** 2
        DK0(I) = UFR ** 2 / CU ** .5 * .5 * (FE / FVI) ** (-.5)
     $           * (DFE * FVI - FE * DFVI) / FVI ** 2
        DT0(I) = 2.5 * UFR ** 2 * TS * FVI ** PR0 / (Z(I)* L * G * .4)
        D2T0(I) = 2.5 * UFR ** 2 * TS / (Z(I)* L * G * .4) *
     $            FVI ** (PR0-1) * (PR0*DFVI - FVI/Z(I))
        D2K0(I) = UFR ** 2 / CU ** .5 * .5 * (FE / FVI) ** (-.5) *
     $            (-.5*DFE**2/FE/FVI + 1.5*FE*DFVI**2/FVI**3 +
     $            D2FE/FVI - (DFE*DFVI+FE*D2FVI)/FVI**2)
        D2U0(I) = 2.5 * UFR * (DFVI * Z(I) - FVI) / Z(I) ** 2
        D2E0(I) = 2.5 * UFR ** 3 / Z(I) ** 4 * (D2FE * Z(I) ** 3 -
     $           (DFE * Z(I) - FE) * 2. * Z(I))

        END DO

        RETURN
        END

C ==========================================================================


        SUBROUTINE CALCULO (PL, NY, NZ, N, ITM, DX, DY, DZ, DENS, BETA,
     $         G, ERRM, ERRORP, RELAJ, PR, SK, SE, CE2, CE1, U, V, W, P,
     $         T, E, K, VI, U0, E0, VI0, DU0, DT0, DK0, DE0, DVI0,
     $         D2U0, D2E0, D2K0, D2T0, K0, U1, V1, W1, T1, E1, K1, U2,
     $         V2, W2, T2, E2, K2, X, A, B, C, D, PP, A1, NP, ND,
     $         DATA$, maxpp, TETA, UT, CASO, DTIEMPO, KT, ET, TIEMPO,
     $         S, MODIFICAR)

C  ----------------------------------------------------------------------
C     Subrutina que resuelve las ecuaciones diferenciales del modelo, y
C  calcula las magnitudes U, V, W, T, P, K, E, VI.
C  ----------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES EXTERNAS

        INTEGER*4 N, ITM, PL, NY, NZ, A1, NP, ND, CASO, S
        INTEGER*4 MODIFICAR(S,7)

        REAL*8 DX, DY, DZ, DENS, BETA, G, ERRM, ERRORP, RELAJ, PR, SK,
     $         SE, CE2, CE1, maxpp, TETA, DTIEMPO, TIEMPO

        REAL*8 U(NY, NZ), V(NY, NZ), W(NY, NZ),
     $         P(NY, NZ), T(NY, NZ), E(NY, NZ),
     $         K(NY, NZ), VI(NY, NZ), UT(NY, NZ, NP+1),
     $         KT(NY, NZ, NP+1), ET(NY, NZ, NP+1)

        REAL*8 U0(NZ), K0(NZ), E0(NZ), VI0(NZ), DU0(NZ), DT0(NZ),
     $         DK0(NZ), DE0(NZ), DVI0(NZ), D2U0(NZ), D2E0(NZ), D2K0(NZ),
     $         D2T0(NZ)

        CHARACTER*40 DATA$

C   VARIABLES INTERNAS

        INTEGER*4 I, J, IT, J1, I1, pantalla, pantallalimite, NTUR, Z1,
     $            Y, kk, IJ

        REAL*8 DD, AP, CP, Er, EP, PP1

        REAL*8 U1(NY, NZ), V1(NY, NZ), W1(NY, NZ),
     $         T1(NY, NZ), E1(NY, NZ), K1(NY, NZ)

        REAL*8 U2(NY, NZ), V2(NY, NZ), W2(NY, NZ),
     $         T2(NY, NZ), E2(NY, NZ), K2(NY, NZ)

        REAL*8 U3(NY, NZ), V3(NY, NZ), W3(NY, NZ),
     $         T3(NY, NZ), E3(NY, NZ), K3(NY, NZ)

        REAL*8 PP(NY, NZ), X(NY+1, NZ+1)

        REAL*8 A(NY - 1), B(NY - 1), C(NY - 1), D(NY - 1)

        pantallalimite=200
        pantalla=1

C   COMIENZO

        DO J = 1 , NY
        DO I = 1 , NZ
             VI(J, I) = .03 * (K0(I) + K(J, I)) ** 2 / (E0(I) + E(J, I))
     $                  - .03 * K0(I) ** 2 / E0(I)
!-----------------------------------
! CAMBIO DE LA VISCOSIDAD TURBULENTA Y ALMACENAMIENTO DE LAS VARIABLES INICIALES
! DENTRO DEL ROTOR EXPANDIDO DURANTE UNA DISTANCIA
! TDWA*(DIAMETRO DE TURBINA) PARA IMPEDIR EL DESARROLLO LIBRE DE LA ESTELA
           DO NTUR=1,S,1
              IF ((MODIFICAR(NTUR,1).LE.N).AND.
     $            (N.LE.MODIFICAR(NTUR,2))) THEN
                 IF ((MODIFICAR(NTUR,5).LE.J).AND.
     $               (J.LE.MODIFICAR(NTUR,6))) THEN
                    Z1 = MODIFICAR(NTUR,7) - 1
                    DO Y = MODIFICAR(NTUR,5),MODIFICAR(NTUR,6),1
                       IF (Y.LE.(MODIFICAR(NTUR,3)-MODIFICAR(NTUR,7)))
     $                          Z1 = Z1 + 1
                       IF (Y.GT.(MODIFICAR(NTUR,3)+MODIFICAR(NTUR,7)))
     $                          Z1 = Z1 - 1
                       kk=MODIFICAR(NTUR,4)-z1
                       IF (kk.LE.0) kk=1
                       IF ((J.EQ.Y).AND.(kk.LE.I).AND.
     $         (I.LE.(MODIFICAR(NTUR,4) + Z1))) then
!                          VI(J, I) = 0.
                          U3(J, I) = U(J, I)
                          V3(J, I) = V(J, I)
                          W3(J, I) = W(J, I)
                          T3(J, I) = T(J, I)
                          K3(J, I) = K(J, I)
                          E3(J, I) = E(J, I)
                       endif
                    END DO
                 ENDIF
              ENDIF
           ENDDO
!-----------------------------------
        END DO
        END DO

        DO J = 1 , NY
        DO I = 1 , NZ
          U1(J, I) = U(J, I)
          V1(J, I) = V(J, I)
          W1(J, I) = W(J, I)
          T1(J, I) = T(J, I)
          K1(J, I) = K(J, I)
          E1(J, I) = E(J, I)
        END DO
        END DO


C   ****EQUATION 1. U COMPONENT****

        DO J = 2 , (NY - 1)
        DO I = 2 , (NZ - 1)
            A(I - 1) = -(VI0(I) + VI(J, I) + (VI(J, I + 1) -
     $                 VI(J, I - 1)) / 4.) / DZ ** 2 - DVI0(I) / DZ / 2
     $                 + W1(J, I) / DZ / 2
            B(I - 1) = 2. * ((U0(I) + U1(J, I)) / DX + (VI0(I) +
     $                 VI(J, I)) / DZ ** 2)
            IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                     B(I - 1) = B(I - 1) + 1. / DTIEMPO
            C(I - 1) = -(VI0(I) + VI(J, I) - (VI(J, I + 1) -
     $                 VI(J, I - 1)) / 4.) / DZ ** 2 + DVI0(I) / DZ / 2
     $                 - W1(J, I) / DZ / 2
            D(I - 1) = 2. * (U0(I) + U1(J, I)) * U1(J, I) / DX - W1(J,I)
     $                 * DU0(I) + (VI0(I) + VI(J, I)) * (U1(J + 1, I)- 2
     $                 * U1(J, I) + U1(J - 1, I)) / DY ** 2 + (VI(J+1,I)
     $                 - VI(J - 1, I)) / DY / 2. * (U1(J + 1, I) -
     $                 U1(J - 1, I)) / DY / 2. + VI(J,I) * D2U0(I) +
     $                 (VI(J, I + 1) - VI(J, I - 1)) * DU0(I) / DZ / 2.
     $                 - V1(J, I) * (U1(J + 1, I) - U1(J - 1, I))/DY/2.
            IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                      D(I - 1) = D(I - 1) +  UT(J,I,N+1) / DTIEMPO
        END DO
           CALL TRIDI (NY, NZ, J, U, A, B, C, D)
        END DO

        DO J = 1 , NY
        DO I = 1 , NZ
           U2(J, I) = U(J, I)
        END DO
        END DO

        DO I = 2 , (NZ - 1)
        DO J = 2 , (NY - 1)
            A(J - 1) = -(VI0(I) + VI(J, I)) / DY ** 2 - (VI(J + 1, I)
     $                 - VI(J - 1, I)) / DY ** 2 / 4. + V1(J, I) /DY/2.
            B(J - 1) = 2. * ((U0(I) + U1(J, I)) / DX + (VI0(I) +
     $                 VI(J, I)) / DY ** 2)
            IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                      B(J - 1) = B(J - 1) + 1. / DTIEMPO
            C(J - 1) = -(VI0(I) + VI(J, I)) / DY ** 2 + (VI(J + 1, I)
     $                 - VI(J - 1, I)) / DY ** 2 / 4. - V1(J, I)/DY/2.
            D(J - 1) = 2. * (U0(I) + U1(J, I)) * U2(J, I)/DX-W1(J, I)
     $                 * DU0(I) + (VI0(I) + VI(J, I)) * (U2(J,I+1) -2
     $                 * U2(J, I) + U2(J, I - 1)) / DZ ** 2 + DVI0(I)
     $                 * (U2(J, I + 1) - U2(J, I - 1)) / DZ / 2. +
     $                 (VI(J, I + 1) - VI(J, I - 1)) * (U2(J, I + 1)
     $                 - U2(J,I-1)) / DZ ** 2 / 4.+VI(J, I) * D2U0(I)
     $                 + DU0(I) * (VI(J, I + 1) - VI(J, I - 1)) /DZ/2.
     $                 - W1(J, I) * (U2(J, I + 1) - U2(J, I - 1))/DZ/2.
            IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                      D(J - 1) = D(J - 1) +  UT(J,I,N+1) / DTIEMPO
        END DO
           CALL TRIDJ (NY, NZ, I, U, A, B, C, D)
        END DO

        IF (CASO.EQ.1) THEN
           DO J = 1 , NY
           DO I = 1 , NZ
              UT(J, I, N+1) = U(J, I)
              V(J, I) = U0(I) * TAN (TETA * 3.141592654 / 180)
           END DO
           END DO
        END IF

        IF (PL.EQ.0) GO TO 357

C   ****EQUATION 2. V COMPONENT****

88      DO J = 2 , (NY - 1)
        DO I = 2 , (NZ - 1)
            A(I - 1) = -(VI0(I) + VI(J, I) + (VI(J, I + 1) - VI(J,
     $                I - 1)) / 4.) / DZ ** 2 - DVI0(I) / DZ / 2. +
     $                W1(J, I) / DZ / 2.
            B(I - 1) = 2. * ((U0(I) + U1(J, I)) / DX + (VI0(I) +
     $                 VI(J, I)) / DZ ** 2)
            C(I - 1) = -(VI0(I) + VI(J, I) - (VI(J, I + 1) - VI(J,
     $                 I - 1)) / 4.) / DZ ** 2 + DVI0(I) / DZ / 2.
     $                 - W1(J, I) / DZ / 2.
            D(I - 1) = 2. * (U0(I) + U1(J, I)) * V1(J, I) / DX -
     $                 (P(J + 1, I) - P(J - 1, I)) / DENS / DY / 2 +
     $                 (VI0(I) + VI(J, I)) * (V1(J + 1, I) - 2. *
     $                 V1(J, I) + V1(J - 1, I)) / DY ** 2 + (W1(J + 1
     $                 , I) - W1(J - 1, I)) * ((VI(J, I + 1) - VI(J,
     $                 I - 1)) / DZ / 2. + DVI0(I)) / DY / 2. +
     $                 (VI(J + 1, I) - VI(J - 1, I)) * (V1(J + 1, I)
     $                 - V1(J - 1, I)) / DY ** 2 / 2. - (K1(J + 1, I)
     $                 - K1(J - 1, I)) / DY / 3. - V1(J, I) * (V1(J + 1
     $                 , I) - V1(J - 1, I)) / DY / 2.
        END DO
            CALL TRIDI(NY, NZ, J, V, A, B, C, D)
        END DO

        DO J = 1 , NY
        DO I = 1 , NZ
          V2(J, I) = V(J, I)
        END DO
        END DO

        DO I = 2 , (NZ - 1)
        DO J = 2 , (NY - 1)
            A(J - 1) = -(VI0(I) + VI(J, I) + (VI(J + 1, I) -
     $                 VI(J - 1, I)) / 2.) / DY ** 2 + V1(J, I) / DY /2.
            B(J - 1) = 2. * ((U0(I) + U1(J, I))/DX + (VI0(I) + VI(J,I))
     $                 / DY ** 2)
            C(J - 1) = -(VI0(I) + VI(J, I)-(VI(J + 1, I) - VI(J - 1, I))
     $                 / 2.) / DY ** 2 - V1(J, I) / DY / 2.
            D(J - 1) = 2. * (U0(I) + U1(J, I)) * V2(J, I) / DX-(P(J+1,I)
     $                 - P(J - 1, I))/DENS/DY/2. + (VI0(I) + VI(J, I))
     $                 * (V2(J, I + 1) - 2. * V2(J, I) + V2(J, I - 1))
     $                 / DZ ** 2 + (V2(J, I + 1) - V2(J, I - 1)) *
     $                 ((VI(J, I + 1) - VI(J, I - 1)) / DZ /2.+DVI0(I))
     $                 / DZ / 2. + (W1(J + 1, I) - W1(J - 1, I)) *
     $                 ((VI(J, I + 1) - VI(J, I - 1)) / DZ / 2.+DVI0(I))
     $                 / DY / 2. - (K1(J + 1, I) - K1(J - 1, I)) / DY/3.
     $                 - W1(J, I) * (V2(J, I + 1) - V2(J,I-1))/DZ/2.
        END DO
            CALL TRIDJ(NY, NZ, I, V, A, B, C, D)
        END DO



C   ****EQUATION 3.  W COMPONENT****

        DO J = 2 , (NY - 1)
        DO I = 2 , (NZ - 1)
           A(I - 1) = -(VI0(I) + VI(J, I) + (VI(J,I+1)-VI(J,I-1)) / 2.)
     $                / DZ ** 2 - DVI0(I) / DZ + W1(J, I) / DZ / 2.
           B(I - 1) = 2. * ((U0(I) + U1(J, I)) / DX+(VI0(I)+VI(J,I))
     $                / DZ ** 2)
           C(I - 1) = -(VI0(I) + VI(J, I) - (VI(J, I + 1)- VI(J,I-1))
     $                / 2.) / DZ ** 2 + DVI0(I) / DZ - W1(J, I)/DZ/2.
           D(I - 1) = 2. *(U0(I)+U1(J, I))*W1(J, I) / DX - (P(J, I + 1)
     $                -P(J,I-1)) / DENS / DZ / 2. + (VI0(I) + VI(J,I))
     $                * (W1(J+1,I) - 2. * W1(J,I) + W1(J-1,I)) / DY ** 2
     $                + (VI(J+1,I) - VI(J-1,I))*(W1(J+1,I) - W1(J-1,I))
     $                / DY ** 2 / 4. + (VI(J+1,I) - VI(J-1,I)) *
     $                (V1(J,I+1) - V1(J,I-1)) /DZ /DY /4. - (K1(J,I+1)
     $                - K1(J,I-1)) /DZ /3. + BETA *G * T1(J,I) - V1(J,I)
     $                * (W1(J + 1, I) - W1(J - 1, I)) / DY / 2.
        END DO
           CALL TRIDI(NY, NZ, J, W, A, B, C, D)
        END DO

        DO J = 1 , NY
        DO I = 1 , NZ
            W2(J, I) = W(J, I)
        END DO
        END DO

        DO I = 2 , (NZ - 1)
        DO J = 2 , (NY - 1)
             A(J - 1) = -(VI0(I) + VI(J,I) + (VI(J+1,I)- VI(J-1,I)) /4.)
     $                  / DY ** 2 + V1(J, I) / DY / 2.
             B(J - 1) = 2. * ((U0(I) + U1(J,I)) /DX + (VI0(I) + VI(J,I))
     $                  / DY ** 2)
             C(J - 1) = -(VI0(I) + VI(J,I)- (VI(J+1,I) - VI(J-1,I)) /4.)
     $                  / DY ** 2 - V1(J, I) / DY / 2.
             D(J - 1) = 2. * (U0(I) + U1(J,I)) * W2(J,I) /DX - (P(J,I+1)
     $                  - P(J,I-1)) /DENS /DZ /2. + (VI0(I) + VI(J,I)) *
     $                  (W2(J,I+1) - 2. * W2(J,I) + W2(J,I-1)) / DZ ** 2
     $                  +(W2(J,I+1)-W2(J,I-1))* ((VI(J,I+1) - VI(J,I-1))
     $                  /DZ /2. + DVI0(I)) /DZ + (VI(J+1,I) - VI(J-1,I))
     $                  * (V1(J,I+1) - V1(J,I-1)) /DY /DZ /4. + BETA* G*
     $                  T1(J,I)-(K1(J,I+1)-K1(J,I-1)) /DZ /3.-W1(J,I) *
     $                  (W2(J, I + 1) - W2(J, I - 1)) / DZ / 2.
        END DO
            CALL TRIDJ(NY, NZ, I, W, A, B, C, D)
        END DO




C  ****** EQUATION 4.  CONTINUITY. PRESSURE CORRECTION ******


        IT = 0
        DD = DY ** 2 + DZ ** 2
        AP = DY ** 2 / DD / 2.
        CP = DZ ** 2 / DD / 2.

180     IT = IT + 1

        Er = 0.

        DO J = 2 , (NY - 1)
        DO I = 2 , (NZ - 1)
             EP = -2. *DENS *(U0(I) + U(J,I)) *AP *CP *DD /DX *((U(J,I)
     $            - U1(J, I)) /DX +(V(J+1,I) - V(J-1,I)) /DY /2. +
     $            (W(J, I + 1) - W(J, I - 1)) / DZ / 2.)
             PP1 = PP(J, I)
             PP(J, I) = AP *(PP(J,I+1) + PP(J,I-1)) + CP * (PP(J + 1, I)
     $                  + PP(J - 1, I)) + EP
             Er = Er + ABS(PP(J, I) - PP1)
        END DO
        END DO

c original       WRITE (*,10000) IT, Er
c original10000   FORMAT (1X, I3, ' ', F7.3)

        maxpp=0.
        DO J = 2 , (NY - 1)
        DO I = 2 , (NZ - 1)
           maxpp=max(maxpp,ABS(PP(J, I)))
        END DO
        END DO

        IF (Er.LT.ERRM) THEN
           if (pantalla.eq.1) then
!              WRITE (*,10001) IT, Er, maxpp, N, NP-N, A1, ND, DATA$
!10001         FORMAT (' Converge in',I4, ' iterations. Errm:',F8.2,
!     $            '(Pa). Errorp: ',F8.2,'(Pa). S-', I3,' Remain:', I3,
!     $            ' Case:', I2,'/', I2, ' (', A40, ')')
              pantalla=pantallalimite
           else
              pantalla=pantalla-1
           endif
        ELSE
           IF (IT.LT.ITM) THEN
             GO TO 180
           ELSE
             if (pantalla.eq.1) then
!               WRITE (*,10002) IT, Er, maxpp, N, NP-N, A1, ND, DATA$
!10002          FORMAT (' Do not converge in',I4, ' iterations. Errm:',
     $                F8.2,'(Pa). Errorp: ',F8.2,'(Pa). S-', I3,
     $            ' Remain:', I3, ' Case:', I2, '/', I2, ' (', A40, ')')
                pantalla=pantallalimite
             else
                pantalla=pantalla-1
             endif
          END IF
        END IF

        DO J = 2 , (NY - 1)
        DO I = 2 , (NZ - 1)
           IF (ABS(PP(J, I)).GT.ERRORP) THEN
              DO J1 = 1 , NY
              DO I1 = 1 , NZ
                 P(J1, I1) = P(J1, I1) + RELAJ * PP(J1, I1)
              END DO
              END DO
              GO TO 88
           END IF
        END DO
        END DO

        DO J = 1 , NY
        DO I = 1 , NZ
           P(J, I) = P(J, I) + PP(J, I)
        END DO
        END DO



C   ****EQUATION 5.  ENERGY****

        DO J = 2 , (NY - 1)
        DO I = 2 , (NZ - 1)
             A(I - 1) = -(VI0(I)+VI(J,I) + (VI(J,I+1) - VI(J,I-1)) / 4.)
     $                  /PR /DZ**2 - DVI0(I) /PR /DZ /2. + W1(J,I)/DZ/2.
             B(I - 1) = 2. *((U0(I) + U1(J,I)) /DX +(VI0(I) + VI(J,I))
     $                  / PR / DZ ** 2)
             C(I - 1) = -(VI0(I) + VI(J,I) - (VI(J,I+1) - VI(J,I-1))/4.)
     $                  /PR /DZ ** 2 + DVI0(I) / PR / DZ / 2. - W1(J, I)
     $                  / DZ / 2.
             D(I - 1) = 2. *(U0(I) + U1(J,I)) * T1(J,I) /DX - W1(J, I)
     $                  * DT0(I) + (VI0(I) + VI(J,I)) /PR * (T1(J+1,I)
     $                  - 2. * T1(J,I) + T1(J-1,I)) /DY**2 + (VI(J+1,I)
     $                  - VI(J - 1, I)) / PR / DY / 2. * (T1(J + 1, I)
     $                  - T1(J-1,I)) /DY /2. + VI(J,I) /PR * D2T0(I)
     $                  + (VI(J,I+1) - VI(J,I-1)) /PR *DT0(I) /DZ /2.
     $                  - V1(J,I) * (T1(J + 1, I)-T1(J-1,I)) / DY / 2.
     $
        END DO
            CALL TRIDI(NY, NZ, J, T, A, B, C, D)
        END DO

        DO J = 1 , NY
        DO I = 1 , NZ
           T2(J, I) = T(J, I)
        END DO
        END DO

        DO I = 2 , (NZ - 1)
        DO J = 2 , (NY - 1)
             A(J - 1) = -(VI0(I) + VI(J,I)) /PR /DY**2 - (VI(J+1,I)
     $                  - VI(J-1,I)) /PR /DY ** 2 /4. + V1(J,I) /DY /2.
             B(J - 1) = 2. * ((U0(I) + U1(J,I)) /DX +(VI0(I) + VI(J,I))
     $                  / PR / DY ** 2)
             C(J - 1) = -(VI0(I) + VI(J,I)) /PR /DY ** 2 +(VI(J+1,I)
     $                  - VI(J-1,I)) /PR /DY ** 2 /4. - V1(J,I) /DY /2.
             D(J - 1) = 2. *(U0(I) + U1(J,I)) * T2(J,I) / DX - W1(J,I)
     $                  * DT0(I) + (VI0(I) + VI(J,I)) /PR *(T2(J,I+1)
     $                  - 2. *T2(J,I) + T2(J,I-1)) /DZ ** 2 + (T2(J,I+1)
     $                  - T2(J,I-1)) * DVI0(I) /PR /DZ /2. +(VI(J, I+1)
     $                  - VI(J,I-1))/PR * (T2(J,I+1)-T2(J,I-1)) /DZ ** 2
     $                  / 4. + VI(J, I) / PR * D2T0(I) + DT0(I) *
     $                  (VI(J,I+1) - VI(J, I - 1)) / PR /DZ /2.
     $                  - W1(J,I) * (T2(J,I+1) - T2(J,I-1)) / DZ / 2.

        END DO
            CALL TRIDJ(NY, NZ, I, T,A , B, C, D)
        END DO



C    ***** EQUATION 6.   TURBULENT KINETIC ENERGY ******

357     DO J = 2 , (NY - 1)
        DO I = 2 , (NZ - 1)
             A(I - 1) = -(VI0(I) + VI(J, I)) / SK / DZ ** 2 - DVI0(I)
     $                  /SK /DZ /2. - (VI(J,I+1) - VI(J,I-1)) /SK /DZ
     $                  ** 2 / 4. + W1(J, I) / DZ / 2.
             B(I - 1) = 2. *(U0(I)+U1(J,I)) /DX +2.*(VI0(I) + VI(J,I))
     $                 /SK /DZ**2 +(E1(J, I) + E0(I))/(K1(J,I)+K0(I))
             IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                       B(I - 1) = B(I - 1) + 1. / DTIEMPO
             C(I - 1) = -(VI0(I) + VI(J, I)) / SK / DZ ** 2 + DVI0(I)
     $                  / SK / DZ / 2. + (VI(J, I + 1) - VI(J, I - 1))
     $                  / SK / DZ ** 2 / 4. - W1(J, I) / DZ / 2.
             D(I - 1) = 2. * (U0(I) + U1(J,I)) / DX * K1(J,I) - W1(J,I)
     $                  * DK0(I) + (VI0(I) + VI(J, I)) /SK * (K1(J+1,I)
     $                  - 2. * K1(J, I) + K1(J - 1, I)) / DY ** 2 +
     $                  (VI(J + 1, I) - VI(J - 1, I)) / SK * (K1(J+1,I)
     $                  - K1(J-1,I)) /DY ** 2 /4. + DK0(I) * (VI(J,I+1)
     $                  - VI(J, I - 1)) / SK / DZ / 2. + VI(J, I) / SK
     $                  * D2K0(I) + (VI0(I) + VI(J, I)) *
     $                  (((U1(J, I + 1) - U1(J, I - 1)) / DZ / 2.) ** 2
     $                  + ((U1(J + 1, I) - U1(J - 1, I)) / DY/2.)**2 +
     $                  2. * DU0(I) * (U1(J,I+1)-U1(J,I-1)) / DZ / 2.)
     $                  - BETA * G / PR * ((VI0(I) + VI(J,I))
     $                  * (T1(J,I+1) - T1(J,I-1)) / DZ / 2. + VI(J, I)
     $                  * DT0(I)) + VI(J, I) * DU0(I) ** 2 - K0(I) *
     $                  (E1(J, I) + E0(I)) / (K1(J, I) + K0(I)) + E0(I)
     $                  - V1(J,I) * (K1(J+1,I) - K1(J-1,I)) / DY / 2.
            IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                      D(I - 1) = D(I - 1) +  KT(J,I,N+1) / DTIEMPO
        END DO
            CALL TRIDI(NY, NZ, J, K, A, B, C, D)
        END DO

        DO J = 1 , NY
        DO I = 1 , NZ
           K2(J, I) = K(J, I)
        END DO
        END DO

        DO I = 2 , (NZ - 1)
        DO J = 2 , (NY - 1)
             A(J - 1) = -(VI0(I) + VI(J,I)) /SK /DY ** 2 - (VI(J+1,I)
     $                  - VI(J - 1, I)) / SK / DY ** 2 / 4 + V1(J, I)
     $                  / DY / 2.
             B(J - 1) = 2. * (U0(I) + U1(J, I)) / DX + 2. * (VI0(I) +
     $                  VI(J, I)) / SK / DY ** 2 + (E1(J, I) + E0(I))
     $                  / (K1(J, I) + K0(I))
             IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                       B(J - 1) = B(J - 1) + 1. / DTIEMPO
             C(J - 1) = -(VI0(I) + VI(J, I)) /SK /DY ** 2 + (VI(J+1,I)
     $                  - VI(J - 1, I)) / SK / DY ** 2 / 4 - V1(J, I)
     $                  / DY / 2.
             D(J - 1) = 2. * (U0(I) + U1(J,I)) / DX *K2(J,I) - W1(J,I)
     $                  * DK0(I) + (VI0(I) + VI(J, I)) * (K2(J, I + 1)
     $                  - 2. * K2(J, I) + K2(J, I - 1)) / SK / DZ ** 2
     $                  + (VI(J, I + 1) - VI(J, I - 1)) * (K2(J, I + 1)
     $                  - K2(J, I - 1)) / SK / DZ ** 2 / 4. + DVI0(I) *
     $                  (K2(J, I + 1) - K2(J, I - 1)) / SK / DZ / 2. +
     $                  DK0(I) / SK * (VI(J, I + 1) - VI(J, I - 1)) / DZ
     $                  / 2. + VI(J, I) * D2K0(I) / SK + (VI0(I) +
     $                  VI(J,I)) * (((U1(J,I+1)- U1(J,I-1))/DZ/2.) ** 2
     $                  + ((U1(J + 1,I) - U1(J - 1, I)) / DY / 2.) ** 2
     $                  + 2. * DU0(I) * (U1(J,I+1) - U1(J,I-1))/DZ / 2.)
     $                   - BETA * G / PR * ((VI0(I) + VI(J, I)) *
     $                  (T1(J, I + 1) - T1(J, I - 1)) / DZ / 2. +
     $                  VI(J, I) * DT0(I)) + VI(J, I) * DU0(I) ** 2
     $                  - K0(I) * (E1(J,I)+E0(I))
     $                  /(K2(J,I)+ K0(I))+ E0(I) - W1(J,I) * (K2(J,I+1)
     $                  - K2(J, I - 1)) / DZ / 2.
            IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                      D(J - 1) = D(J - 1) +  KT(J,I,N+1) / DTIEMPO
        END DO
            CALL TRIDJ(NY, NZ, I, K, A, B, C, D)
        END DO

        IF (CASO.EQ.1) THEN
           DO J = 1 , NY
           DO I = 1 , NZ
              KT(J, I, N+1) = K(J, I)
           END DO
           END DO
        END IF


C   ****EQUATION 7. DISSIPATION RATE OF K (Epsilon)****

        DO J = 2 , (NY - 1)
        DO I = 2 , (NZ - 1)
             A(I - 1) = -(VI0(I) + VI(J, I)) / SE / DZ ** 2 - DVI0(I)
     $                  / SE / DZ / 2. - (VI(J, I + 1) - VI(J, I - 1))
     $                  / SE / DZ ** 2 / 4. + W1(J, I) / DZ / 2.
             B(I - 1) = 2. * (U0(I) + U1(J, I)) / DX + 2. * (VI0(I) +
     $                  VI(J,I)) /SE /DZ ** 2 +CE2 *(E1(J,I) + E0(I))
     $                  / (K1(J, I) + K0(I))
             IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                       B(I - 1) = B(I - 1) + 1. / DTIEMPO
             C(I - 1) = -(VI0(I) + VI(J, I)) / SE / DZ ** 2 + DVI0(I)
     $                  / SE / DZ / 2. + (VI(J, I + 1) - VI(J, I - 1))
     $                  / SE / DZ ** 2 / 4. - W1(J, I) / DZ / 2.
             D(I - 1) = 2. *(U0(I) + U1(J,I)) * E1(J,I) / DX - W1(J,I)
     $                  * DE0(I) + (VI0(I) + VI(J,I)) /SE * (E1(J+1,I)
     $                  - 2. *E1(J,I) +E1(J-1,I)) /DY ** 2 +(VI(J+1,I)
     $                  - VI(J-1, I)) * (E1(J + 1, I) - E1(J - 1, I))
     $                  /SE /DY ** 2 /4. +DE0(I) /SE *(VI(J,I+1) -
     $                  VI(J,I-1)) /DZ /2. + VI(J,I) / SE *D2E0(I)+CE1
     $                  *(E1(J, I)+E0(I)) /(K1(J, I) +K0(I)) * (VI0(I)
     $                  + VI(J, I)) * (((U1(J, I + 1) - U1(J, I - 1))
     $                  / DZ / 2.) ** 2 + ((U1(J + 1, I) - U1(J-1, I))
     $                  / DY / 2.) ** 2 + 2. * DU0(I) * (U1(J, I + 1)
     $                  - U1(J, I - 1)) / DZ / 2. + DU0(I) ** 2 - .2 *
     $                  BETA * G * (DT0(I) + (T1(J,I+1) - T1(J,I-1))
     $                  / DZ / 2.) / PR) - CE1 * E0(I) / K0(I) * VI0(I)
     $                  * (DU0(I) ** 2 - .2 * BETA * G * DT0(I)) - CE2
     $                  * (E1(J, I) + E0(I)) / (K1(J,I) + K0(I)) * E0(I)
     $                  + CE2*E0(I)**2 /K0(I) -V1(J,I) *(E1(J+1,I) -
     $                  E1(J - 1, I)) / DY / 2.
            IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                      D(I - 1) = D(I - 1) +  ET(J,I,N+1) / DTIEMPO
        END DO
            CALL TRIDI(NY, NZ, J, E, A, B, C, D)
        END DO

        DO J = 1 , NY
        DO I = 1 , NZ
             E2(J, I) = E(J, I)
        END DO
        END DO

        DO I = 2 , (NZ - 1)
        DO J = 2 , (NY - 1)
             A(J - 1) = -(VI0(I) + VI(J,I)) /SE /DY ** 2 - (VI(J+1,I)
     $                  - VI(J - 1, I)) / SE / DY ** 2 / 4 + V1(J, I)
     $                  / DY / 2.
             B(J - 1) = 2.*(U0(I) + U1(J,I)) /DX +2. *(VI0(I) + VI(J,I))
     $                  /SE /DY**2 + CE2 * (E1(J,I) + E0(I)) / (K1(J, I)
     $                  + K0(I))
             IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                       B(J - 1) = B(J - 1) + 1. / DTIEMPO
             C(J - 1) = -(VI0(I) + VI(J, I)) / SE / DY ** 2 + (VI(J+1,I)
     $                  - VI(J - 1, I)) / SE / DY ** 2 / 4 - V1(J, I) /
     $                  DY / 2.
             D(J - 1) = 2. *(U0(I) + U1(J,I)) * E2(J, I) / DX - W1(J, I)
     $                  * DE0(I) + (VI0(I) + VI(J,I)) / SE * (E2(J,I+1)
     $                  - 2. * E2(J,I) + E2(J,I-1)) /DZ ** 2 +(VI(J,I+1)
     $                  - VI(J, I - 1)) * (E2(J, I + 1) - E2(J, I - 1))
     $                  / SE / DZ ** 2 / 4. + DE0(I) / SE * (VI(J, I+1)
     $                  - VI(J,I-1)) / DZ / 2. + VI(J,I) / SE * D2E0(I)
     $                  + CE1 * (E2(J, I) + E0(I)) / (K1(J, I) + K0(I))
     $                  * (VI0(I) + VI(J,I))*(((U1(J,I+1)-U1(J,I-1)) /
     $                  DZ / 2.) ** 2 + ((U1(J + 1, I) - U1(J - 1, I))
     $                  / DY / 2.) ** 2 + 2. * DU0(I) * (U1(J, I + 1) -
     $                  U1(J, I - 1)) / DZ / 2. + DU0(I) ** 2 -.2*BETA
     $                  * G * (DT0(I) + (T1(J, I + 1) - T1(J, I - 1)) /
     $                  DZ / 2.) / PR) - CE1 * E0(I) / K0(I) * VI0(I) *
     $                  (DU0(I) ** 2 - .2 * BETA * G * DT0(I)) - CE2 *
     $                  (E2(J, I) + E0(I)) / (K1(J,I) + K0(I)) * E0(I)
     $                  + CE2 * E0(I) ** 2/K0(I) - W1(J,I) * (E2(J,I+1)
     $                  - E2(J, I - 1)) / DZ / 2.
            IF ((CASO.EQ.1).AND.(TIEMPO.GT.0.))
     $                      D(J - 1) = D(J - 1) +  ET(J,I,N+1) / DTIEMPO
        END DO
            CALL TRIDJ(NY, NZ, I, E, A, B, C, D)
        END DO

        IF (CASO.EQ.1) THEN
           DO J = 1 , NY
           DO I = 1 , NZ
              ET(J, I, N+1) = E(J, I)
           END DO
           END DO
        END IF

        DO J = 2 , (NY - 1)
        DO I = 2 , (NZ - 1)
             VI(J, I) = .03 * (K0(I) + K(J, I)) ** 2 / (E0(I) + E(J, I))
     $                   - .03 * K0(I) ** 2 / E0(I)
!-----------------------------------
! CAMBIO DE LA VISCOSIDAD TURBULENTA E IMPOSICI�N DE LAS VARIABLES INICIALES DEL C�LCULO
! COMO VARIABLES DE SALIDA DENTRO DEL ROTOR EXPANDIDO DURANTE UNA DISTANCIA
! TDWA*(DIAMETRO DE TURBINA) PARA IMPEDIR EL DESARROLLO LIBRE DE LA ESTELA
           DO NTUR=1,S,1
              IF ((MODIFICAR(NTUR,1).LE.N).AND.
     $            (N.LE.MODIFICAR(NTUR,2))) THEN
                 IF ((MODIFICAR(NTUR,5).LE.J).AND.
     $               (J.LE.MODIFICAR(NTUR,6))) THEN
                    Z1 = MODIFICAR(NTUR,7) - 1
                    DO Y = MODIFICAR(NTUR,5),MODIFICAR(NTUR,6),1
                       IF (Y.LE.(MODIFICAR(NTUR,3)-MODIFICAR(NTUR,7)))
     $                          Z1 = Z1 + 1
                       IF (Y.GT.(MODIFICAR(NTUR,3)+MODIFICAR(NTUR,7)))
     $                          Z1 = Z1 - 1
                       kk=MODIFICAR(NTUR,4)-z1
                       IF (kk.LE.0) kk=1
                       IF ((J.EQ.Y).AND.(kk.LE.I).AND.
     $         (I.LE.(MODIFICAR(NTUR,4) + Z1))) then
!                          VI(J, I) = 0.
                          U(J, I) = U3(J, I)
                          V(J, I) = V3(J, I)
                          W(J, I) = W3(J, I)
                          T(J, I) = T3(J, I)
                          K(J, I) = K3(J, I)
                          E(J, I) = E3(J, I)
                          U1(J, I) = U3(J, I)
                          V1(J, I) = V3(J, I)
                          W1(J, I) = W3(J, I)
                          T1(J, I) = T3(J, I)
                          K1(J, I) = K3(J, I)
                          E1(J, I) = E3(J, I)
                       endif
                    END DO
                 ENDIF
              ENDIF
           ENDDO
!-----------------------------------
        END DO
        END DO

        RETURN
        END

C ==========================================================================


        SUBROUTINE CALCULOKMWM(CONTROLKMWM,NSKMWM,SECKMWM,
     $                N,NY,NZ,DELTAV0,U,DELTAVM,DTIEMPO,DELTAK)

C -----------------------------------------------------------
C  CALCULO DE ANALISIS DEL Kinematic Model For Wake Meandering
C -----------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES EXTERNAS

        INTEGER*4 N, NY, NZ, CONTROLKMWM, NSKMWM

        INTEGER*4 SECKMWM(NSKMWM)

        REAL*8 DTIEMPO,U(NY, NZ),DELTAV0(NSKMWM),
     +         DELTAVM(NSKMWM,NY,NZ),DELTAK(NSKMWM,NY,NZ)

C   VARIABLES INTERNAS

        INTEGER*4 IJ, J, I

C   COMIENZO

           DO IJ=1,NSKMWM,1
              IF (SECKMWM(IJ).EQ.N) THEN
                 IF (CONTROLKMWM.EQ.1) THEN
                    DO J = 1 , NY
                       DO I = 1 , NZ
                          DELTAV0(IJ)=MIN(DELTAV0(IJ),U(J, I))
                       END DO
                    END DO
                    ELSEIF (CONTROLKMWM.EQ.2) THEN
                       DO J = 1 , NY
                          DO I = 1 , NZ
                             DELTAVM(IJ ,J, I)=
     $                            DELTAVM(IJ, J, I)+U(J, I)*DTIEMPO
                          END DO
                       END DO
                       ELSEIF (CONTROLKMWM.EQ.3) THEN
                          DO J = 1 , NY
                             DO I = 1 , NZ
                                DELTAK(IJ ,J, I)=
     $  DELTAK(IJ, J, I)+((U(J, I)-DELTAVM(IJ, J, I))**2.)*DTIEMPO*0.5
                             END DO
                          END DO
                          ELSE
!                             WRITE(*,*)'NO SE EMPLEA ESTE C�LCULO'
                             STOP
                 ENDIF
              ENDIF
           ENDDO

        RETURN
        END

C ==========================================================================

        SUBROUTINE CAMBIOEJES (S, NTS, X2, Y2, X1, Y1, XTS, YTS, XTS1,
     $                         YTS1, AB1)

C  -----------------------------------------------------------
C  Subrutina para el cambio de los ejes arbitrarios iniciales
C  a otros ejes (X2,Y2,Z2), en los que el eje X2 tiene la direc-
C  cion del viento
C  -----------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I

        REAL*8 AB

C   VARIABLES EXTERNAS

        INTEGER*4 S, NTS

        REAL*8 AB1

        REAL*8 X2(S), Y2(S), X1(S), Y1(S)

        REAL*8 XTS(NTS), YTS(NTS), XTS1(NTS), YTS1(NTS)


C   COMIENZO

        AB = AB1 * 3.141592654 / 180

        DO I = 1 , S
           X2(I) = X1(I) * COS(AB) - Y1(I) * SIN(AB)
           Y2(I) = Y1(I) * COS(AB) + X1(I) * SIN(AB)
        END DO

        DO I = 1 , NTS
           XTS1(I) = XTS(I) * COS(AB) - YTS(I) * SIN(AB)
           YTS1(I) = YTS(I) * COS(AB) + XTS(I) * SIN(AB)
        END DO

        RETURN
        END

C ==========================================================================


        SUBROUTINE CONTORNO (S, I, NY, NZ, DIAM, CT, DY, HY, HZ, VM,
     $                 VMAX, U0, U, VM1, T, TK, K, AK, E, K0, E0,
     $                 JINF, JSUP, GR)

C  ------------------------------------------------------------------
C    Subrutina que calcula la velocidad incidente en el rotor para
C  la estimacion de la potencia que produce la aeroturbina.
C  ------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 GRID, HY1, HY2, Z1, NUM, Y, Z ,kk

        REAL*8 DIAMEXP

C   VARIABLES EXTERNAS

        INTEGER*4 I, HY, HZ, S, NY, NZ, JINF, JSUP, GR

        REAL*8 CT, DY, AK, TK, INCU, VM1

        REAL*8 DIAM (S), VM(S), VMAX(S), U0(NZ), U(NY, NZ), T(NY,NZ),
     $         K(NY, NZ), E(NY, NZ), K0(NZ), E0(NZ)


C   CALCULO DEL DEFECTO DE VELOCIDAD EN EL ROTOR EXPANDIDO DE INTERACCION

        INCU = VM1 * (1 - SQRT(1 - CT))

C   ESTIMACION DEL DIAMETRO EXPANDIDO.

        DIAMEXP = DIAM(I) * SQRT((1 + SQRT(1 - CT)) / 2 / SQRT(1 - CT))
        GRID = INT(DIAMEXP / DY / 2)

C   CAMBIO DE LAS CONDICIONES DE CONTORNO EN EL ROTOR EXPANDIDO.
C   CALCULO DE LA VELOCIDAD INCIDENTE EN EL ROTOR EXPANDIDO SEGUN EL
C   FLUJO BASICO PARA LA DETERMINACION DEL RENDIMIENTO.

        IF ((INT(GRID / 2)).EQ.(GRID / 2)) THEN
           GR = (GRID - 1) / 2
        ELSE
           GR = GRID / 2
        END IF

        HY1 = HY - GRID
        HY2 = HY + GRID
        Z1 = GR - 1
        NUM = 0
        VM(I) = 0.
        VMAX(I) = 0.
        JINF=HY1
        JSUP=HY2

        DO Y = HY1 , HY2

           IF (Y.LE.(HY - GR)) THEN
             Z1 = Z1 + 1
           END IF
           IF (Y.GT.(HY + GR)) THEN
             Z1 = Z1 - 1
           END IF
		 kk=hz-z1
		 IF (kk.LE.0) kk=1

           DO Z =kk , (HZ + Z1)
CFerm�n     DO Z = (HZ - Z1) , (HZ + Z1)
             VM(I) = VM(I) + (U0(Z) + U(Y, Z)) ** 3
             VMAX(I) = VMAX(I) + (U0(Z)) ** 3
             NUM = NUM + 1
             U(Y, Z) = U(Y, Z) - INCU           ! NUEVAS CONDICIONES EN
             T(Y, Z) = T(Y, Z) + TK             ! EL ROTOR EXPANDIDO
             K(Y, Z) = K(Y, Z) + AK
             E(Y, Z) = K(Y, Z) / K0(Z) * E0(Z)
           END DO
        END DO

        VM(I) = (VM(I) / NUM) ** (1. / 3.)

        VM(I) = VM1 ! CON ESTA SENTENCIA OBLIGO A QUE SE EMPLEE LA VELOCIDAD EN EL CENTRO DEL ROTOR PARA EL C�LCULO DE LAS POSICIONES DE LA TURBINA

        VMAX(I) = (VMAX(I) / NUM) ** (1. / 3.)

        RETURN
        END


C ==========================================================================


        SUBROUTINE CORRECVEL (ND, NT, A1, DI, NM, AB1, NM1, DIR1, VH2,
     $                        VH, DIR, OROG$, VH1, DIREC, MAST, NTMAX,
     $                        DIRMAX)

C  -----------------------------------------------------------------------
C       Subrutina que corrige la velocidad medida en el mastil segun la
C   influencia de la orografia en el emplazamiento del mismo; a esta
C   velocidad corregida estan referidos el resto de los porcentajes
C   de influencia de la orografia en los emplazamientos de las turbinas.
C  -----------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I, D

        REAL*8 MAST1

C   VARIABLES EXTERNAS

        INTEGER*4 NT, A1, NM1, DIR1, ND, NTMAX, DIRMAX

        INTEGER*4 NM(DIR1), DIR (ND)

        REAL*8 AB1, VH1

        REAL*8 DI(ND), VH2(ND), VH(ND)

        REAL*8 DIREC(DIRMAX), MAST (NTMAX, DIRMAX)

        CHARACTER*40 OROG$

C     CASO EN QUE NO SE TENGA EN CUENTA LA OROGRAFIA

        VH2(A1) = VH(A1)

C     CASO EN QUE SI SE TENGA EN CUENTA LA OROGRAFIA

        IF (OROG$.NE.'no') THEN

C          AVERIGUAMOS SI LA DIRECCION ACTUAL COINCIDE CON ALGUNA DE LAS
C          ANTERIORES

           IF (NT.NE.1) THEN
             IF (A1.GT.1) THEN
                DO I = 1 , (A1 - 1)
                  IF (INT(DI(A1)).EQ.INT(DI(I))) THEN
                     NM(A1) = NM(I)
                     WRITE (*,10000) NM(A1)
10000                FORMAT (' NUMBER OF THE MAST THAT MEASURES',
     $                         ' WIND INCIDENT SPEED = ',I3)
                     GO TO 2
                  END IF
                END DO
             END IF
3            WRITE (*,10001) AB1
10001        FORMAT (/,1X,'WIND DIRECTION = ',F6.2,'�')
             WRITE (*,*) 'NUMBER OF THE MAST THAT MEASURES WIND ',
     $                   'INCIDENT SPEED ? '
             READ (*,'(I3)') NM(A1)
             IF (NM(A1).GT.NT) THEN
                WRITE (*,10003) NT
10003           FORMAT (/, 1X,'THERE ARE ONLY ', I3, ' MASTS.')
                GO TO 3
             END IF
2            NM1 = NM(A1)
           ELSE
             NM1 = 1
           END IF

C        CALCULO DE LA VELOCIDAD INCIDENTE EN EL PARQUE (CORRECCION DE LA
C        MEDIDA EN EL MASTIL METEOROLOGICO).

           DO D = 1 , (DIR1 - 1)
             IF ((DI(A1).GE.DIREC(D)).AND.(DI(A1).LT.DIREC(D + 1))) THEN
               MAST1 = MAST(NM1, D + 1) + (DI(A1) - DIREC(D + 1)) *
     $                (MAST(NM1, D) - MAST(NM1, D + 1)) / (DIREC(D) -
     $                DIREC(D + 1))
               VH2(A1) = 100 * VH(A1) / (100 + MAST1)
               DIR(A1) = D
             END IF
           END DO

           IF (DI(A1).GE.DIREC(DIR1)) THEN
             MAST1 = MAST(NM1, 1) + (DI(A1)-DIREC(1))*(MAST(NM1,DIR1)
     $         - MAST(NM1, 1)) / (DIREC(DIR1) - DIREC(1))
             VH2(A1) = 100 * VH(A1) / (100 + MAST1)
             DIR(A1) = DIR1
           END IF
        END IF

        VH1 = VH2(A1)

        RETURN
        END

C ==========================================================================


        SUBROUTINE ENCABEC (ND, OROG$, DATA$, NM1, ALT, VH, AB1,
     $                     FIC, A1, TIEMPO)

C  --------------------------------------------------------------------
C     Subrutina que a�ade al fichero de potencias y rendimientos el
C  encabezamiento que separa cada caso de analisis. Tambien lo indica
C  por pantalla.
C  --------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS


C   VARIABLES EXTERNAS

        INTEGER*4 A1, NM1, ND, FIC

        REAL*8 ALT, AB1, TIEMPO

        REAL*8 VH(ND)

        CHARACTER*40 OROG$, DATA$

C   COMIENZO DE SUBRUTINA

        IF (A1.EQ.1) THEN

          IF (FIC.EQ.7) THEN
            WRITE (FIC,10000)
10000       FORMAT (' U velocity component',//)
          END IF
          IF (FIC.EQ.9) THEN
            WRITE (FIC,10001)
10001       FORMAT (' W velocity component',//)
          END IF
          IF (FIC.EQ.8) THEN
            WRITE (FIC,10002)
10002       FORMAT (' V velocity component',//)
          END IF
          IF (FIC.EQ.11) THEN
            WRITE (FIC,10003)
10003       FORMAT (' Temperature T',//)
          END IF
          IF (FIC.EQ.10) THEN
            WRITE (FIC,10004)
10004       FORMAT (' Pressure P',//)
          END IF
          IF (FIC.EQ.12) THEN
            WRITE (FIC,10005)
10005       FORMAT (' Turbulent kinetic energy K',//)
          END IF
          IF (FIC.EQ.13) THEN
            WRITE (FIC,10006)
10006       FORMAT (' Epsilon E',//)
          END IF
          IF (FIC.EQ.14) THEN
            WRITE (FIC,10007)
10007       FORMAT (' Turbulent vicosity VI',//)
          END IF
          IF (FIC.EQ.17) THEN
            WRITE (FIC,10020)
10020       FORMAT (' Module of Velocity',//)
          END IF
          IF (FIC.EQ.2) THEN
            WRITE (FIC, 9998)
9998        FORMAT (' BASIC FLOW ',//)
          END IF
          IF (FIC.EQ.15) THEN
            WRITE (FIC, 9999)
9999        FORMAT (' TURBULENT SPECTRA ',//,' Turbulence spectra'
     $              ,' in J/kg for different frecuencies',//)
          END IF
          IF (FIC.NE.6) THEN
            WRITE (FIC, 10008) DATA$
10008       FORMAT (1X, 'PARK DATA FILE: ',A,/)
          END IF

        END IF



        WRITE (FIC, 10009)
10009   FORMAT (/, ' -----------------------------',
     $           '--------------------------- ')
        IF (OROG$.NE.'no') THEN
           WRITE (FIC, 10010) NM1
10010      FORMAT ('     MAST THAT MEASURES WIND INCIDENT',
     $              ' SPEED = ', I3)
           WRITE (FIC, 10011) ALT
10011      FORMAT ('     MAST HEIGHT = ', F6.2,' m')
           WRITE (FIC, 10012) VH(A1)
10012      FORMAT ('     WIND SPEED MEASURED IN MAST = ',
     $              F6.2,' m/s')
           WRITE (FIC, 10013) AB1
10013      FORMAT ('     WIND DIRECTION MEASURED IN MAST = ',
     $              F6.2,'�')
           WRITE (FIC, 10050) TIEMPO
10050      FORMAT ('     TIME = ', F8.2,'s')
           IF (FIC.EQ.90) THEN
              WRITE(FIC, 23456)
23456         FORMAT ('     SE SUPONE QUE LA TURBINA DE MAYOR ',
     $'�NDICE ES EL M�STIL DE MEDICI�N, POR TANTO, PARA EL C�LCULO ',
     $'DE LAS CONDICIONES AMBIENTALES EN DICHA POSICI�N SE USA LA ',
     $'ALTURA DEL M�STIL METEOROL�GICO',/)
           ENDIF
        ELSE
           WRITE (FIC, 10014) ALT
10014      FORMAT ('     MAST HEIGHT = ', F6.2, ' m')
           WRITE (FIC, 10015) VH(A1)
10015      FORMAT ('     UNPERTURBED WIND SPEED MEASURED IN',
     $           ' MAST = ', F6.2, ' m/s')
           WRITE (FIC, 10016) AB1
10016      FORMAT ('     WIND DIRECTION MEASURED IN MAST = ',
     $              F6.2, '�')
           WRITE (FIC, 10051) TIEMPO
10051      FORMAT ('     TIME = ', F8.2,'s')
           IF (FIC.EQ.90) THEN
              WRITE(FIC, 12345)
12345         FORMAT ('     SE SUPONE QUE LA TURBINA DE MAYOR ',
     $'�NDICE ES EL M�STIL DE MEDICI�N, POR TANTO, PARA EL C�LCULO ',
     $'DE LAS CONDICIONES AMBIENTALES EN DICHA POSICI�N SE USA LA ',
     $'ALTURA DEL M�STIL METEOROL�GICO')
           ENDIF
        END IF
        WRITE (FIC, 10017)
10017   FORMAT (' -----------------------------------------',
     $           '--------------- ',/)

        IF (FIC.EQ.15) THEN
           WRITE (FIC, 10018)
10018      FORMAT (/,'   X ,  Y ,  Z (m);    n = 0.001     0.01',
     $            '     0.1       1        10    (1/s)')
           WRITE (FIC, 10019)
10019      FORMAT ('  ---- ---- ----       --- -----    ------',
     $            '   -----    -----    -----',/)
        END IF

        RETURN
        END
C ==========================================================================


        FUNCTION nSu (AA,BB)

C ----------------------------------------------------------------
C     Funcion que calcula el espectro de la turbulencia segun la
C energia cinetica turbulenta y la frecuencia adimensionalizada.
C ----------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

C   VARIABLES EXTERNAS

        REAL*8 BB, AA, nSu

C   COMIENZO DE FUNCION

        nSu = 19.2*AA*BB / (1. + 33.*AA)**(5./3.)

        RETURN
        END

C ==========================================================================


        SUBROUTINE ESPECTRO (XTS, YTS, ZTS, XTS1, YTS1, DX, DY, DZ, NTS,
     $                     DISTAN, U, K, E, U0, K0, E0, NY, NZ, YMINMA,
     $                     S, MINX, FNC)

C  ----------------------------------------------------------------
C       Subrutina que calcula el espectro de la turbulencia en los
C   puntos especificados en el fichero de definicion del parque.
C  ----------------------------------------------------------------


C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 HZ1, HY1, I, J, D1, D2, D3

        REAL*8 KK, EE, ZZ, FF(5), TT(5)

C   VARIABLES EXTERNAS

        INTEGER*4 NTS, NY, NZ, S

        REAL*8 DX, DY, DZ, DISTAN, YMINMA, MINX

        REAL*8 XTS(NTS), YTS(NTS), ZTS(NTS), XTS1(NTS), YTS1(NTS),
     $         U0(NZ), K0(NZ), E0(NZ), U(NY,NZ), K(NY,NZ), E(NY,NZ)

        REAL*8 FNC

C   COMIENZO SUBRUTINA

c original        D1 = INT((MINX + DISTAN - DX) * 10000)
c original        D2 = INT((MINX + DISTAN) * 10000)
        D1 = INT((DISTAN - DX) * 10000)
        D2 = INT((DISTAN) * 10000)

        DO I = 1 , NTS
c original           D3 = INT(XTS1(I) * 10000)
           D3 = INT((XTS1(I) - MINX) * 10000)
           IF ((D3.GT.D1).AND.(D3.LE.D2)) THEN
              HZ1 = INT(ZTS(I) / DZ + .5) + 1
              HY1 = INT((YTS1(I) - YMINMA) / DY + .5) + 1
              KK = K0(HZ1)+K(HY1, HZ1)
              EE = E0(HZ1)+E(HY1, HZ1)
              ZZ = 2.5*(KK/5.47)**1.5 / EE
              DO J = 1 , 5
                FF(J) = 10.**(J-4) * ZZ / U0(HZ1)
                TT(J) = FNC(FF(J),KK)
              END DO
              WRITE (15, 10002) XTS(I), YTS(I), ZTS(I),(TT(J), J= 1, 5)
10002         FORMAT (1X, 3F6.2,6X,5(F7.4,2X),/)
           END IF
        END DO

        RETURN
        END

C ==========================================================================


       SUBROUTINE FICHOROG (OROG$, S, DIR1, NT, VO, DIREC, MAST, VMED,
     $                      DIRMAX, NTMAX)

C -----------------------------------------------------------------
C    Subrutina que lee el fichero de orografia, formado a partir de
C los datos del codigo WASP. Procesa los datos correspondientes.
C -----------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES EXTERNAS

        INTEGER*4 S, NT, DIR1, DIRMAX, NTMAX

        REAL*8 VO(S, DIRMAX), DIREC(DIRMAX), MAST (NTMAX, DIRMAX)

        CHARACTER*40 OROG$

C   VARIABLES INTERNAS

        INTEGER*4 M, T, D

        REAL*8 VMED(DIRMAX)

        CHARACTER*80 B$


C    LECTURA DEL FICHERO DE OROGRAFIA OROG$

        OPEN (16, FILE = OROG$)
          READ (16, '(I3)') DIR1
          READ (16,'(I3)') NT
          READ (16,'(A)') B$
          DO D = 1 , DIR1
             READ (16, *) DIREC(D), (VO(T,D), T=1,S),
     $                        (MAST(M,D), M=1,NT)
          END DO
        CLOSE (16)


C    CALCULO DEL VALOR MEDIO DE LAS CORRECCIONES PARA CADA DIRECCION

        DO D = 1 , DIR1
           VMED(D) = 0.
           DO M = 1 , NT
              VMED(D) = MAST(M, D) + VMED(D)
           END DO
           DO T = 1 , S
              VMED(D) = VMED(D) + VO(T, D)
           END DO
           VMED(D) = VMED(D) / (S + NT)
        END DO

C    REFERIMOS LA TABLA DE OROGRAFIA AL VALOR MEDIO PARA CADA DIRECCION

        DO D = 1 , DIR1
           DO T = 1 , S
              VO(T, D) = VO(T, D) - VMED(D)
           END DO
           DO M = 1 , NT
              MAST(M, D) = MAST(M, D) - VMED(D)
           END DO
        END DO

        RETURN
        END

C ==========================================================================

        SUBROUTINE GUARDARKMWM(NY,NZ,NSKMWM,SECKMWM,
     +                         DY,DZ,U,DELTAV0,DELTAVM,DELTAK)

C  -----------------------------------------------------------------------
C  ALMACENAMIENTO DE LOS RESULTADOS DE KMWM
C  -----------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 IJ,J,I

C   VARIABLES EXTERNAS

        INTEGER*4 NY,NZ,NSKMWM
        INTEGER*4 SECKMWM(NSKMWM)

        REAL*8  DY,DZ
	  REAL*8  U(NY,NZ),DELTAV0(NSKMWM)
	  REAL*8  DELTAVM(NSKMWM,NY,NZ),DELTAK(NSKMWM,NY,NZ)

C   COMIENZO SUBRUTINA

           OPEN(100,FILE='KMWM.dat')
           DO IJ=1,NSKMWM,1
              WRITE(100,FMT='(2A)')' COMPROBACI�N DE KINEMATIC MODEL ',
     +                           'FOR WAKE MEANDERING'
              WRITE(100,FMT='(2A,I3.3)')' CALCULO DE DELTAK/DELTAV0^2 ',
     +                                'PARA LA SECCI�N: ',SECKMWM(IJ)

              WRITE (100,*)
              WRITE (100, 20004)' DELTAV0: ',DELTAV0(IJ)
20004         FORMAT (A,F9.4)

              WRITE (100,*)
              WRITE (100,*)'DISTANCIA EN Y'
              WRITE (100, 20001)
     +             (J*DY,J=(NY - 1),2,-1)
20001         FORMAT (3000F9.4)

              WRITE (100,*)
              WRITE (100,*)'DELTAVM'
              DO I = 1, NZ
                 WRITE (100, 20002)
     +             I*DZ,(DELTAVM(IJ ,J, I),J=(NY - 1),2,-1)
20002            FORMAT (3000F9.4)
              END DO

              WRITE (100,*)
              WRITE (100,*)'DELTAK'
              DO I = 1, NZ
                 WRITE (100, 20003)
     +             I*DZ,(DELTAK(IJ ,J, I),J=(NY - 1),2,-1)
20003            FORMAT (3000F9.4)
              END DO

              WRITE (100,*)
              WRITE (100,*)'DELTAK/DELTAV0^2'
              DO I = 1, NZ
                 WRITE (100, 20000)
     +         I*DZ,(DELTAK(IJ ,J, I)/(DELTAV0(IJ))**2.,J=(NY - 1),2,-1)
20000            FORMAT (3000F9.4)
              END DO
              WRITE (100,*)
              WRITE (100,FMT='(2A)')" **********************************
     $********************************************"
           ENDDO
           ENDFILE(100)
           CLOSE(100)

        RETURN
        END


C ==========================================================================


        SUBROUTINE IMPBSFLOW (NZ, Z, U0, K0, E0, VI0)

C  -----------------------------------------------------------------------
C      Subrutina que almacena en el fichero '.BSF' las magnitudes de velo-
C  cidad V, energia cinetica turbulenta K, disipacion de la energia cine-
C  tica turbulenta E y viscosidad turbulenta VI, correspondientes al flujo
C  basico y calculadas en el programa principal.
C  -----------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I

C   VARIABLES EXTERNAS

        INTEGER*4 NZ

        REAL*8  Z(NZ), U0(NZ), K0(NZ), E0(NZ), VI0(NZ)

C   COMIENZO SUBRUTINA


        WRITE (2, 10008)
10008   FORMAT ( /, 1X, '                 U(m/s)',
     $          '       K(m2/s2)     E(m2/s3)   VI(m2/s)')
        DO I = NZ , 1 , -1
           WRITE (2, 10009) I,Z(I),U0(I),K0(I),E0(I),VI0(I)
10009      FORMAT (1X, I3,'  ',F7.2,' m   ',F8.3,'     ',F8.3,
     $                 '     ', F8.3, '   ', F8.3)
        END DO

        RETURN
        END
C ==========================================================================



        SUBROUTINE INCREMEN (S, DIAM1, DIAM, DY, DZ)

C  ----------------------------------------------------------------
C     Subrutina que determina los incrementos DY, DZ de la malla de
C  calculo del parque.
C  ----------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I, FLAG1

        REAL*8 INTER

C   VARIABLES EXTERNAS

        INTEGER*4 S

        REAL*8 DY, DZ

        REAL*8 DIAM(S), DIAM1(S)


C   COMIENZO SUBRUTINA

C   DETERMINACION DEL DIAMETRO MINIMO

        DO I = 1 , S
          DIAM1(I) = DIAM(I)
        END DO

C     ORDENACION DE LOS DIAMETROS DE MENOR A MAYOR

        FLAG1 = 1
        DO WHILE (FLAG1.EQ.1)
          FLAG1 = 0
          DO I = 1 , (S - 1)
            IF (DIAM1(I).GT.DIAM1(I + 1)) THEN
              INTER = DIAM1(I)
              DIAM1(I) = DIAM1(I + 1)
              DIAM1(I + 1) = INTER
              FLAG1 = 1
            END IF
          END DO
        END DO

C   INCREMENTOS DY y DZ RESPECTO AL DIAMETRO MINIMO

        DZ = DIAM1(1) / 6.
        DY = DZ

C       En teoria, asignando a cada nudo su cuadrado correspondiente en la
C       malla, para representar la superficie del rotor correctamente habria
C       que dividir el diametro en siete partes; pero para tener en cuenta
C       el diametro expandido, es mejor dividirlo por seis.


        RETURN
        END

C ==========================================================================


        SUBROUTINE INICMALL (NY, NZ, U, V, W, P, PP, T, K, E)


C  ----------------------------------------------------------------
C     Subrutina que determina los incrementos DY, DZ de la malla de
C  calculo del parque.
C  ----------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I, J

C   VARIABLES EXTERNAS

        INTEGER*4 NY, NZ

        REAL*8 U(NY, NZ), V(NY, NZ), W(NY, NZ), PP(NY, NZ),
     $         P(NY, NZ), T(NY, NZ), E(NY, NZ), K(NY, NZ)

C   COMIENZO SUBRUTINA

        DO J = 1 , NY
        DO I = 1 , NZ

          U(J, I) = 0.      ! VELOCIDAD LONGITUDINAL
          V(J, I) = 0.      ! VELOCIDAD SEGUN EJE Y
          W(J, I) = 0.      ! VELOCIDAD SEGUN EJE Z
          P(J, I) = 0.      ! PRESION
          PP(J, I) = 0.     ! CORRECCION DE LA PRESION
          T(J, I) = 0.      ! TEMPERATURA
          K(J, I) = 0.      ! ENERGIA CINETICA TURBULENTA
          E(J, I) = 0.      ! DISIPACION DE ENERGIA CINETICA TURBULENTA

        END DO
        END DO

        RETURN
        END

C ==========================================================================


        SUBROUTINE INTERAC (S, MINX, DX, DISTAN, X2, FLAG, FLAG2)

C  ------------------------------------------------------------
C       Subrutina que averigua si entre las dos ultimas secciones
C   calculadas se encuentra alguna turbina.
C  ------------------------------------------------------------


C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I, D1, D2, D3

C   VARIABLES EXTERNAS

        INTEGER*4 S, FLAG2

        INTEGER*4 FLAG(S)

        REAL*8 DX, DISTAN, MINX

        REAL*8 X2(S)


C   COMIENZO SUBRUTINA

        DO I = 1 , S
           FLAG(I) = 0
        END DO

c original        D1 = INT((MINX + DISTAN - DX) * 10000)
c original        D2 = INT((MINX + DISTAN) * 10000)
        D1 = INT((DISTAN - DX) * 10000)
        D2 = INT((DISTAN) * 10000)

        FLAG2 = 0
        DO I = 1 , S
c original           D3 = INT(X2(I)*10000)
           D3 = INT((X2(I) - MINX)*10000)
           IF ((D3.GT.D1).AND.(D3.LE.D2)) THEN
              FLAG(I) = 1
              FLAG2 = 1
           END IF
        END DO

        RETURN
        END


C ==========================================================================


        SUBROUTINE LIMMALLA (X2, Y2, X3, Y3, XTS1, S, NTS, DIAM1, DY,
     $             NY, NSL, DX, NP, HT, Z0, DZ, NZ, YMINMA, MINX, NDIAM)

C  -------------------------------------------------------------------
C      Subrutina que calcula los limites de la malla de calculo y el
C  numero de nodos en los ejes X (NP), Y (NY) y Z (NZ).
C  -------------------------------------------------------------------


C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I, FLAG1

        REAL*8 INTER, YMAXMA, XMA, ZMA

C   VARIABLES EXTERNAS

        INTEGER*4 S, NY, NZ, NP, NSL, NTS

        REAL*8 DX, DY, DZ, YMINMA, HT, Z0, MINX, NDIAM

        REAL*8 X2(S), X3(S+NTS), Y2(S), Y3(S), DIAM1(S), XTS1(NTS)


C   COMIENZO SUBRUTINA

        DO I = 1 , S
           Y3(I) = Y2(I)
           X3(I) = X2(I)
        END DO

        DO I = 1 , NTS
           X3(I+S) = XTS1(I)
        END DO

C       ORDENACION DE LAS TURBINAS A LO LARGO DEL EJE Y

        FLAG1 = 1
        DO WHILE (FLAG1.EQ.1)
           FLAG1 = 0
           DO I = 1 , (S - 1)
             IF (Y3(I).GT.Y3(I + 1)) THEN
               INTER = Y3(I)
               Y3(I) = Y3(I + 1)          ! ORDENACION DE MENOR A MAYOR
               Y3(I + 1) = INTER
               FLAG1 = 1
             END IF
           END DO
        END DO

C       ORDENACION DE LAS TURBINAS A LO LARGO DEL EJE X

        FLAG1 = 1
        DO WHILE (FLAG1.EQ.1)
           FLAG1 = 0
           DO I = 1 , (S + NTS - 1)
             IF (X3(I).GT.X3(I + 1)) THEN
               INTER = X3(I)              ! ORDENACION DE MENOR A MAYOR
               X3(I) = X3(I + 1)
               X3(I + 1) = INTER
               FLAG1 = 1
             END IF
           END DO
        END DO

C   NUMERO DE NODOS EN LA MALLA SEGUN EL EJE Y

!        YMINMA = Y3(1) - 3 * DIAM1(S)     ! SE A�ADEN POR AMBOS EXTREMOS
!        YMAXMA = Y3(S) + 3 * DIAM1(S)     ! TRES DIAMETROS.
        YMINMA = Y3(1) - NDIAM * DIAM1(S)     ! SE A�ADEN POR AMBOS EXTREMOS
        YMAXMA = Y3(S) + NDIAM * DIAM1(S)     ! NDIAM DIAMETROS.

        NY = INT((YMAXMA - YMINMA) / DY) + 1

C   NUMERO DE SECCIONES SEGUN EL EJE X

        IF (NSL.EQ.0) THEN
          NSL = 1
        END IF
        XMA = X3(S+NTS) - X3(1)
        MINX = X3(1)
        NP = INT(XMA / DX) + NSL

C   CALCULO DEL NUMERO DE NODOS EN EL EJE Z

!        ZMA = HT + 2.5 * DIAM1(S)      ! SE A�ADEN DOS DIAMETROS Y MEDIO
        ZMA = HT + NDIAM * DIAM1(S)      ! SE A�ADEN NDIAM DIAMETROS
        NZ = INT((ZMA - Z0) / DZ) + 1  ! POR ENCIMA DE LA ALTURA MEDIA.

        RETURN
        END


C ==========================================================================

        SUBROUTINE MEANDERING(mx,my,mz,ponde,per,R,
     $                        dtiempo,tiempolimite,DATA$,ntur)

C ----------------------------------------------------------------------
C      Subrutina que promedia los distintos casos de meandering.
C ----------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 i,j,kk,cm,paso,salir,casos,salto,continuar,turbina

        real*8 vb,db,vts
        real*8, allocatable, dimension(:) :: matrizaux,vt
        real*8, allocatable, dimension(:,:,:) :: matrizu,matrizk

        CHARACTER*40 ficherou,ficherok,ficheropow
        CHARACTER*1 a1,aux1
        character*2 a2
        character*4 a4
        character*27 a27
        character*39 a39
        character*42 a42
        character*45 a45
        character*46 a46
        character*48 a48
        character*57 a57
        character*58 a58
        character*19 a19
        CHARACTER*10 formato10
        CHARACTER*50 nombre13
        character*60 nombreu,nombrek,nombrepow


C   VARIABLES EXTERNAS

        INTEGER*4 mx,my,mz,ponde,per,R,ntur

        REAL*8 dtiempo,tiempolimite

        CHARACTER*40 DATA$

C COMIENZO DE SUBRUTINA

           open(1,file='INDICE.DAT')
           read(1,*)aux1
           ficherou=DATA$(:R) // '.gru'
           ficherok=DATA$(:R) // '.grk'
           ficheropow=DATA$(:R) // '.pow'
           read(1,*)aux1
           paso=tiempolimite/dtiempo
           salir=0
           casos=0
           do while(salir.eq.0)
              read(1,*,err=8812)aux1
              casos=casos+1
           enddo
 8812      close(1)

           if (ponde.eq.0) then
              salto=casos-paso
              if (casos.ne.(per*paso+1)) then
                 continuar=2
                 do while ((continuar.ne.0).and.(continuar.ne.1))
                    write(*,fmt='(2a)')
     +                  'No se ponderan pasos exactos ',
     +                  'de aire por el parque.'
                    write(*,fmt='(a,f8.2,a,f8.2,a)')
     +                  'El paso real es de ',
     +                  tiempolimite,'s, pero se va a promediar con ',
     +                  paso*dtiempo,'s'
                    write(*,fmt='(a)')'�Quiere continuar? (1=si, 0=no)'
                    read(*,*)continuar
                 enddo
                 if (continuar.eq.0) stop
              endif
              write(*,*)
     +'***********************************************************'
              write(*,*)
     +'SE PROMEDIA �NICAMENTE EL �LTIMO PASO DE AIRE POR EL PARQUE'
              write(*,*)
     +'***********************************************************'
              pause
              elseif (ponde.eq.1) then
                 salto=0
                 paso=casos
           endif

           allocate(vt(1:ntur))
           allocate(matrizaux(1:my),
     +              matrizu(1:mx,1:my,1:mz),
     +              matrizk(1:mx,1:my,1:mz))
           matrizu=0.
           matrizk=0.
	     vt=0.
           write(formato10,fmt='(a1,i4,a5)')'(',my,'f7.2)'

           open(4,file=ficherou)
           open(5,file=ficherok)
           open(1,file='INDICE.DAT')
           read(1,*)aux1
           read(1,*)aux1

           do cm=1,casos,1

              read(1,fmt='(a)')nombre13
              if (cm.le.salto) cycle !�nicamente se elige el �ltimo paso de aire por el parque para la ponderaci�n para asegurar estacionareidad en el proceso

              write(*,fmt='(a,i3.3,a,i3.3,2a)')'Procesando el fichero ',
     +                           cm-salto,' de ',paso,' ',nombre13
              nombreu(1:r+9)=nombre13(1:r+9)
              nombreu(r+10:r+13)='.gru'
              nombrek(1:r+9)=nombre13(1:r+9)
              nombrek(r+10:r+13)='.grk'
              nombrepow(1:r+9)=nombre13(1:r+9)
              nombrepow(r+10:r+13)='.pow'


              open(2,file=nombreu)
              open(3,file=nombrek)

              if (cm.eq.(salto+1)) then
                 read(2,fmt='(a27)')a27
                 write(4,fmt='(a27)')a27
                 read(2,fmt='(a2)')a2
                 write(4,fmt='(a2)')a2
                 read(2,fmt='(a2)')a2
                 write(4,fmt='(a2)')a2
                 read(2,fmt='(a45)')a45
                 write(4,fmt='(a45)')a45
                 read(2,fmt='(a2)')a2
                 write(4,fmt='(a2)')a2
                 read(2,fmt='(a2)')a2
                 write(4,fmt='(a2)')a2
                 read(2,fmt='(a58)')a58
                 write(4,fmt='(a58)')a58
                 read(2,fmt='(a27)')a27
                 write(4,fmt='(a27)')a27
                 read(2,fmt='(a57)')a57
                 write(4,fmt='(a57)')a57
                 read(2,fmt='(a46)')a46
                 write(4,fmt='(a46)')a46
                 read(2,fmt='(a19)')a19
                 write(4,fmt='(a19)')a19
                 read(2,fmt='(a58)')a58
                 write(4,fmt='(a58)')a58
                 write(4,fmt='(a1)')' '

                 read(3,fmt='(a27)')a27
                 write(5,fmt='(a27)')a27
                 read(3,fmt='(a2)')a2
                 write(5,fmt='(a2)')a2
                 read(3,fmt='(a2)')a2
                 write(5,fmt='(a2)')a2
                 read(3,fmt='(a45)')a45
                 write(5,fmt='(a45)')a45
                 read(3,fmt='(a2)')a2
                 write(5,fmt='(a2)')a2
                 read(3,fmt='(a2)')a2
                 write(5,fmt='(a2)')a2
                 read(3,fmt='(a58)')a58
                 write(5,fmt='(a58)')a58
                 read(3,fmt='(a27)')a27
                 write(5,fmt='(a27)')a27
                 read(3,fmt='(a57)')a57
                 write(5,fmt='(a57)')a57
                 read(3,fmt='(a46)')a46
                 write(5,fmt='(a46)')a46
                 read(3,fmt='(a19)')a19
                 write(5,fmt='(a19)')a19
                 read(3,fmt='(a58)')a58
                 write(5,fmt='(a58)')a58
                 write(5,fmt='(a1)')' '
                 else
                    do i=1,12,1
                       read(2,fmt='(a)')a2
                       read(3,fmt='(a)')a2
                    enddo
              endif

              do i=1,mx,1
        	      do kk=1,mz,1
                   matrizaux=0.
                   read(2,fmt=formato10,err=1234)
     +                    (matrizaux(my+1-j),j=1,my)
1234               do j=1,my,1
                      matrizu(i,j,kk)=matrizu(i,j,kk)+matrizaux(j)/paso
                   enddo
                enddo
              enddo

              do i=1,mx,1
                 do kk=1,mz,1
                    matrizaux=0.
                    read(3,fmt=formato10,err=4567)
     +                     (matrizaux(my+1-j),j=1,my)
4567                do j=1,my,1
                       matrizk(i,j,kk)=matrizk(i,j,kk)+matrizaux(j)/paso
                    enddo
                 enddo
              enddo

              close(2)
              close(3)

              open(2,file=nombrepow)
         do i=1,3,1
            read(2,*)a1
         enddo
         read(2,fmt='(a48,f5.2)')a48,vb
         read(2,fmt='(a39,f6.2)')a39,db
         do i=1,2,1
            read(2,*)a1
         enddo
         do i=1,2,1
            read(2,*)
         enddo
         do i=1,ntur,1
            read(2,fmt='(a42,i1.1,a4,f5.2)')a42,turbina,a4,vts
            vt(i)=vt(i)+vts/paso
            read(2,*)
         enddo
              close(2)


           enddo

           close(1)

           write(*,*)'ESCRIBIENDO LOS FICHEROS DE SALIDA'

           do i=1,mx,1
              do kk=1,mz,1
                 write(4,fmt=formato10)(matrizu(i,my+1-j,kk),j=1,my,1)
                 write(5,fmt=formato10)(matrizk(i,my+1-j,kk),j=1,my,1)
              enddo
           enddo

           endfile(4)
           close(4)
           endfile(5)
           close(5)

           open(4,file=ficheropow)
            write(4,fmt='(a1,//)')a1
         do i=1,5,1
            write(4,fmt='(a1)')a1
         enddo
            write(4,fmt='(a1,//)')a1
         do i=1,ntur,1
            write(4,fmt='(a42,i1.1,a4,f5.2,/)')a42,i,a4,vt(i)
         enddo
           endfile(4)
           close(4)

           deallocate(matrizaux,matrizu,matrizk,vt)

        RETURN
        END

C ==========================================================================


        SUBROUTINE OROGPOT (I, A1, DIR, DIR1, AB1, VM, VMAX, VH1, S,
     $                      VO, DIREC, ND, DIRMAX)

C ----------------------------------------------------------------------
C      Subrutina que estima una velocidad Uorog que resumira el efecto
C la orografia en el emplazamiento de la turbina en cuestion. Se supon-
C dra que dicho efecto se superpone linealmente al de las estelas.
C ----------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 D

        REAL*8 Uorog, VO1

C   VARIABLES EXTERNAS

        INTEGER*4 S, A1, DIR1, ND, I, DIRMAX

        INTEGER*4 DIR(ND)

        REAL*8 AB1, VH1

        REAL*8 VM(S), VMAX(S)

        REAL*8 VO(S, DIRMAX), DIREC(DIRMAX)

C COMIENZO DE SUBRUTINA

        D = DIR(A1)

        IF (D.LT.DIR1) THEN
          VO1 = VO(I, D + 1) + (AB1 - DIREC(D + 1)) * (VO(I, D) -
     $         VO(I, D + 1)) / (DIREC(D) - DIREC(D + 1))
        ELSE
          VO1 = VO(I, 1) + (AB1 - DIREC(1)) * (VO(I, D) - VO(I, 1))
     $         / (DIREC(D) - DIREC(1))
        END IF

        Uorog = VO1 * VH1 / 100  !EFECTO DE LA OROGRAFIA EN LA VELOCIDAD

        VM(I) = VM(I) + Uorog
        VMAX(I) = VMAX(I) + Uorog


        RETURN
        END

C ==========================================================================


        SUBROUTINE POTENCIA (S, VM, FLAG, PCU, POT$, FICH, WARD)

C  -------------------------------------------------------------------
C     Subrutina para el calculo de la potencia de cada aeroturbina,
C  segun la velocidad cubica media calculada en el programa principal
C  y la correccion de la orografia Uorog, si es que la orografia se
C  tiene en cuenta en el analisis del parque.
C  -------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I, T, NDAT

        REAL*8 VA, VD, POT, V1, V2, POT1, POT2

        CHARACTER*8 B$
        CHARACTER*12 FORMATO

C   VARIABLES EXTERNAS

        INTEGER*4 S

        INTEGER*4 WARD(S), FICH(S), FLAG(S)

        REAL*8 VM(S), PCU(S)

        CHARACTER*16 POT$(S)

C COMIENZO DE SUBRUTINA

        DO I = 1 , S

           IF (FLAG(I).NE.0) THEN

           IF (POT$(I).EQ.'POWER CURVE') THEN
             T = 1
             rewind(fich(i))
             WRITE(FORMATO,fmt='(a1,I2.2,a9)')'(',WARD(I),'(/),f6.2)'
             READ (FICH(I),FORMATO) VA
             READ (FICH(I),'(F6.2)') VD
             IF ((VM(I).LT.VA).OR.(VM(I).GT.VD)) THEN
                POT = 0.
                GO TO 642
             END IF

C            DETERMINACION DEL VALOR DE LA POTENCIA POR INTERPOLACION
C            EN SU GRAFICA, SEGUN EL TIPO DE TURBINA QUE CORRESPONDA.

             READ (FICH(I), '(I3)') NDAT
             READ (FICH(I), '(A)') B$
             READ (FICH(I), *) V1, POT1
             IF (VM(I).LT.V1) THEN
                 POT = (-POT1) / (VA - V1) * (VM(I) - V1) + POT1
             ELSE
641              T = T + 1
                 IF (T.GT.NDAT) THEN
                   POT = POT1
                 ELSE
                   READ (FICH(I), *) V2, POT2
                   IF ((VM(I).GE.V1).AND.(VM(I).LT.V2)) THEN
                      POT = (POT2-POT1)/(V2-V1)*(VM(I)-V1)+POT1
                   ELSE
                      V1 = V2
                      POT1 = POT2
                      GO TO 641
                   END IF
                 END IF
             END IF
           END IF

C      IMPRESION EN PARQUE$

642        PCU(I) = POT
           IF (POT$(I).EQ.'POWER CURVE') THEN
             WRITE (3, 10000) I, PCU(I), I, VM(I)
10000        FORMAT (/, 1X, 'POWER IN TURBINE ', I3, ' = ', F7.2,
     $           ' kW     V', I3, ' = ', F6.2, ' m/s')
           ELSE
             WRITE (3, 10001) I, I, VM(I)
10001        FORMAT (/, 1X, 'POWER CURVE IN TURBINE ', I2,
     $              ' NOT AVAILABLE   V', I2, ' = ', F6.2, ' m/s')
           END IF

        END IF

        END DO

        RETURN
        END

C ==========================================================================


        SUBROUTINE RENDIMIENTOS (PCU, S, VMAX, TURB$, COMP, fich)

C -------------------------------------------------------------------------
C    Subrutina para el calculo de la potencia total del parque y de los
C rendimientos de las aeroturbinas y del parque; para su determinacion
C consideraremos que la potencia maxima producida por una aeroturbina en
C el caso de analisis en cuestion corresponde a la velocidad cubica media
C en el rotor de la misma segun el flujo basico.
C -------------------------------------------------------------------------


C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I, T, NDAT, A, FLAG

        REAL*8 VA, VD, POT, V1, V2, POT1, POT2, PPARK, P, RT, RP, B

        CHARACTER*16 POT1$

        CHARACTER*80 B$

        logical*4  od

C   VARIABLES EXTERNAS

        INTEGER*4 S
        integer*4 fich(s)

        REAL*8 COMP

        REAL*8 PCU(S), VMAX(S)

        CHARACTER*25 TURB$(S)

C COMIENZO DE SUBRUTINA

C ****** CALCULO DE LA POTENCIA TOTAL DEL PARQUE ******

        P = 0.
        DO I = 1 , S
          P = P + PCU(I)
        END DO

C       IMPRESION EN PARQUE$ (LABEL 3)

        WRITE (3, 10000) P/1000.
10000   FORMAT (/,1X, 'PARK TOTAL POWER = ',F8.1,' MW',//)

        PPARK = 0.
        FLAG = 0

C ****** CALCULO DE RENDIMIENTOS ******

        DO I = 1 , S

           POT = 0.

C          CALCULO DE LA POTENCIA MAXIMA DE LA TURBINA I POR
C          INTERPOLACION DE LA GRAFICA CORRESPONDIENTE.

           inquire(unit=fich(i),opened=od)
           if(od) then
              REWIND FICH(I)
              else
                 write(*,fmt='(a,i2.2,a)')' Error2: La unidad ',
     +                 fich(i),' no est� abierta'
                 stop
           end if

1002       A = 0
           DO WHILE (A.NE.0)
              read(fich(i),*,err=1003) B
           END DO

1003        read(fich(i),'(A16)') POT1$
           IF ((POT1$.NE.'NO POWER CURVE').AND.(POT1$.NE.'POWER CURVE'))
     $     THEN
              GO TO 1002
           END IF

           IF (POT1$.EQ.'NO POWER CURVE') THEN
              WRITE (3, 10001) I
10001         FORMAT (1X, 'POWER CURVE OF TURBINE ', I3, ' NOT ',
     $               'AVAILABLE',/,1X,'CANNOT CALCULATE EFFICIENCY',/)
              FLAG = 1
              GO TO 1000
           ELSE
              T = 1
              read(fich(i),'(f7.2)') va
              read(fich(i),'(f7.2)') vd
              IF (VMAX(I).LT.VA) THEN
                 VMAX(I) = VA
              END IF
              IF (VMAX(I).GT.VD) THEN
                 VMAX(I) = VD
              END IF
              read (fich(i), '(I3)') NDAT
              read (fich(i), '(A)') B$
              read (fich(i), *) V1, POT1
              IF (VMAX(I).LT.V1) THEN
                  POT = (-POT1) / (VA - V1) * (VMAX(I) - V1) + POT1
              ELSE
1001              T = T + 1
                  IF (T.GT.NDAT) THEN
                     POT = POT1
                  ELSE
                     read(fich(i),*) V2, POT2
                     IF ((VMAX(I).GE.V1).AND.(VMAX(I).LT.V2)) THEN
                       POT = (POT2-POT1)/(V2-V1)*(VMAX(I)-V1)+POT1
                     ELSE
                       V1 = V2
                       POT1 = POT2
                       GO TO 1001
                     END IF
                  END IF
              END IF
           END IF


C          ***** CALCULO DEL RENDIMIENTO DE LA TURBINA I *****

           RT = PCU(I) / POT
           WRITE (3, 10003) I, RT
10003      FORMAT (1X, 'EFFICIENCY TURBINE ', I3,' = ', F4.2,/)

C          CALCULO ACUMULATIVO DE LA MAXIMA POTENCIA QUE PUEDE PRODUCIR
C          EL PARQUE

           PPARK = PPARK + POT

1000       continue
c original1000        close (fich(i))

        END DO


C ***** CALCULO DEL RENDIMIENTO DEL PARQUE *****

        IF ((ABS(PPARK).GT.COMP).AND.(FLAG.EQ.0)) THEN
           RP = P / PPARK
           WRITE (3, 10004) RP
10004      FORMAT (1X,'PARK EFFICIENCY = ', F4.2, /)
        END IF

        RETURN
        END

C ==========================================================================


        SUBROUTINE RESULTADOS (N, NY, NZ, Q, FIC)

C ----------------------------------------------------------------------
C   Subrutina de almacenamiento de resultados segun lo especificado en
C  el fichero de definicion del parque.
C ----------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I, J

C   VARIABLES EXTERNAS

        INTEGER*4 N, NY, NZ, FIC

        REAL*8 Q(NY,NZ)

C COMIENZO DE SUBRUTINA


        IF (FIC.EQ.7) THEN
c           WRITE (FIC, 10010) N
10010      FORMAT (1X, 'U COMPONENT OF VELOCITY. SECTION ',I3,/)
        END IF
        IF (FIC.EQ.9) THEN
c           WRITE (FIC, 10011) N
10011      FORMAT (1X, 'W COMPONENT OF VELOCITY. SECTION ',I3,/)
        END IF
        IF (FIC.EQ.8) THEN
c           WRITE (FIC, 10012) N
10012      FORMAT (1X, 'V COMPONENT OF VELOCITY. SECTION ',I3,/)
        END IF
        IF (FIC.EQ.11) THEN
c           WRITE (FIC, 10013) N
10013      FORMAT (1X, 'TEMPERATURE T. SECTION ' ,I3,/)
        END IF
        IF (FIC.EQ.10) THEN
c           WRITE (FIC, 10014) N
10014      FORMAT (1X, 'PRESSURE P. SECTION ' ,I3,/)
        END IF
        IF (FIC.EQ.12) THEN
c           WRITE (FIC, 10015) N
10015      FORMAT (1X,'TURBULENT KINETIC ENERGY K. SECTION ',I3,/)
        END IF
        IF (FIC.EQ.13) THEN
c           WRITE (FIC, 10016) N
10016      FORMAT (1X, 'EPSILON E. SECTION ',I3,/)
        END IF
        IF (FIC.EQ.14) THEN
c           WRITE (FIC, 10117) N
10117      FORMAT (1X, 'TURBULENT VISCOSITY VI. SECTION ' ,I3,/)
        END IF
        IF (FIC.EQ.17) THEN
c           WRITE (FIC, 10118) N
10118      FORMAT (1X, 'MODULE OF VELOCITY. SECTION ',I3,/)
        END IF

        DO I = 1, NZ
           WRITE (FIC, 10000) (Q(J, I), J = (NY - 1) , 2 , -1)
10000      FORMAT (3000F7.2)
        END DO

c        WRITE (FIC, 10018)
10018   FORMAT (/)

        RETURN
        END


C ==========================================================================

        SUBROUTINE TRIDI (NY, NZ, J, X, A, B, C, D)

C ------------------------------------------------
C   Subrutina de calculo interno del metodo k-E.
C ------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 I

C   VARIABLES EXTERNAS

        INTEGER*4 NY, NZ, J

        REAL*8 X(NY,NZ)

        REAL*8 A(NY - 1), B(NY - 1), C(NY - 1), D(NY - 1)

C COMIENZO DE SUBRUTINA

        DO I = 2 , (NZ - 2)
           B(I) = B(I) - C(I) / B(I - 1) * A(I - 1)
           D(I) = D(I) - C(I) / B(I - 1) * D(I - 1)
        END DO
        X(J, NZ - 1) = D(NZ - 2) / B(NZ - 2)
        DO I = (NZ - 3) , 1 , -1
           X(J, I + 1) = (D(I) - A(I) * X(J, I + 2)) / B(I)
        END DO

        RETURN
        END

C ==========================================================================


        SUBROUTINE TRIDJ (NY, NZ, I, X, A, B, C, D)

C ------------------------------------------------
C   Subrutina de calculo interno del metodo k-E.
C ------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 J

C   VARIABLES EXTERNAS

        INTEGER*4 NY, NZ, I

        REAL*8 X(NY,NZ)

        REAL*8 A(NY - 1), B(NY - 1), C(NY - 1), D(NY - 1)

C COMIENZO DE SUBRUTINA

        DO J = 2 , (NY - 2)
           B(J) = B(J) - C(J) / B(J - 1) * A(J - 1)
           D(J) = D(J) - C(J) / B(J - 1) * D(J - 1)
        END DO
        X(NY - 1, I) = D(NY - 2) / B(NY - 2)
        DO J = (NY - 3) , 1 , -1
           X(J + 1, I) = (D(J) - A(J) * X(J + 2, I)) / B(J)
        END DO

        RETURN
        END

C ==========================================================================


        SUBROUTINE TURBINA (S, I, FICH, AK, TK, VM1, CT, POT$, WARD,
     $                     CTMAX, DIAMCAMBIADO)

C  ----------------------------------------------------------------------
C     Subrutina que lee los datos del fichero de turbina correspondiente
C  y que interpola el coeficiente de empuje Ct.
C  ----------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 T, NDAT, A

        REAL*8 V1, V2, CT1, CT2, B

        CHARACTER*16 CT$, A$

        CHARACTER*80 B$

C     VARIABLES EXTERNAS

        INTEGER*4 I, S

        INTEGER*4 WARD(S), FICH(S)

        REAL*8 AK, TK, VM1, CT, CTMAX, DIAMCAMBIADO

        CHARACTER*16 POT$(S)

         real*8  diam(s)


C COMIENZO DE SUBRUTINA
        rewind(fich(i))
        read(fich(i), '(f6.2)') diam(i)
        READ (FICH(I),'(F6.2)') AK
        READ (FICH(I),'(F6.2)') TK
        READ (FICH(I),'(A)') CT$
        IF (CT$.EQ.'NO CT CURVE') THEN
           WARD(I) = 5
           WRITE (*,10000)
10000      FORMAT (/,1X, 'Ct CURVE NOT AVAILABLE')
           WRITE (*,10001) I, VM1
10001      FORMAT (1X, 'INCIDENT SPEED AT HUB HEIGHT',
     $                ' IN TURBINE ', I3, ' : ',F6.2, ' m/s')
           WRITE (*,*) ' Ct VALUE? '
           READ (*,'(F6.4)') CT
        ELSE
           T = 1                       !DETERMINACION DEL COEFICIENTE DE
           READ (FICH(I), '(I3)') NDAT !PENETRACION Ct POR INTERPOLACION
           READ (FICH(I),*) B$         !EN SU GRAFICA.
           READ (FICH(I), *) V1, CT1
           WARD(I) = NDAT + 7
           IF (VM1.LT.V1) THEN
             WRITE (*,10003)
10003        FORMAT (/,1X,'INCIDENT SPEED OUT OF RANGE IN Ct',
     $             ' CURVE.')
             WRITE (*,10004) I, VM1, V1
10004        FORMAT (1X, 'SPEED AT HUB HEIGHT IN TURBINE ', I3,
     $              ' = ', F6.2, ' m/s  <', F6.2,' m/s (LOW LIMIT)',/)
             WRITE (*,*) ' Ct VALUE? '
             CT = 0.
!             READ (*,'(F6.4)') CT
             WRITE(*,'(A, F6.4)')'CT =', CT
           ELSE
640          T = T + 1
             IF (T.GT.NDAT) THEN
                WRITE (*,10003)
                WRITE (*,10006) I, VM1, V1
10006           FORMAT (1X, 'SPEED AT HUB HEIGHT IN TURBINE ', I3,
     $              ' = ', F6.2, ' m/s  >', F6.2,' m/s (TOP LIMIT)',/)
                WRITE (*,*) ' Ct VALUE? '
                CT = 0.
!                READ (*,'(F6.4)') CT
                WRITE(*,'(A, F6.4)')'CT =', CT
             ELSE
                READ (FICH(I),*) V2, CT2
                IF ((VM1.GE.V1).AND.(VM1.LT.V2)) THEN
                   CT = (CT2 - CT1) / (V2 - V1) * (VM1 - V1) + CT1
                   IF (CT.GT.CTMAX) THEN
                      WRITE (91, *)'==================================='
                      WRITE (91,10007) I, VM1, CT, DIAMCAMBIADO
10007                 FORMAT (1X, 'TURBINE ', I3 , ', V = ', F5.2 ,
     $                        ' m/s, CT = ' , F6.4 ,', DIAMETER = ' ,
     $                        F6.2, ' m ',/)
                      WRITE (91,10008) CT, CTMAX
10008                 FORMAT (1X, 'CT = ', F6.4 , ' GREATER THAN CTMAX',
     $                  ' = ', F6.4 , ' THEREFORE PROPERTIES CHANGE',/)

                      DIAMCAMBIADO = SQRT(CT/CTMAX) * DIAMCAMBIADO
                      CT = CTMAX

                      WRITE (91,10009) I, VM1, CT, DIAMCAMBIADO
10009                 FORMAT (1X, 'TURBINE ', I3 , ', V = ', F5.2 ,
     $                        ' m/s, CT = ' , F6.4 ,', DIAMETER = ' ,
     $                        F6.2, ' m ',/)
                      WRITE (91, *)'==================================='
                   ELSE
                      WRITE (91,10010) I, VM1, CT, DIAMCAMBIADO
10010                 FORMAT (1X, 'TURBINE ', I3 , ', V = ', F5.2 ,
     $                        ' m/s, CT = ' , F6.4 ,', DIAMETER = ' ,
     $                        F6.2, ' m ',/)
                   ENDIF
                ELSE
                   V1 = V2
                   CT1 = CT2
                   GO TO 640
                END IF
             END IF
           END IF
        END IF

        IF ((CT$.EQ.'POWER CURVE').OR.(CT$.EQ.'NO POWER CURVE')) THEN
           POT$(I) = CT$
        ELSE
641        A = 0
           DO WHILE (A.NE.0)
              READ (FICH(I), * , ERR=642) B
           END DO

642        READ (FICH(I),'(A)') A$

           IF ((A$.NE.'NO POWER CURVE').AND.(A$.NE.'POWER CURVE'))
     $     THEN
              GO TO 641
           END IF

           POT$(I) = A$
        END IF

        RETURN
        END

C ==========================================================================


        SUBROUTINE VCEN (VM1, U0, U, HY, HZ, NY, NZ)

C  -----------------------------------------------------------------------
C     Subrutina que calcula la velocidad incidente en el centro del rotor.
C  -----------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

C   VARIABLES EXTERNAS

        INTEGER*4 HY, HZ, NY, NZ

        REAL*8 VM1, U0(NZ), U(NY, NZ)

C   COMIENZO DE SUBRUTINA

C   DETERMINACION DE LA VELOCIDAD EN EL CENTRO DEL ROTOR

        VM1 = U0(HZ) + U(HY, HZ)

        RETURN
        END

C ==========================================================================


        SUBROUTINE KCEN (K, E, HY, HZ, NY, NZ, X1, Y1, I2, S,
     $                   Z, U0, K0, E0, HM, Z2)

C  -----------------------------------------------------------------------
C     Subrutina que almacena la K y E incidente en el centro del rotor.
C  -----------------------------------------------------------------------

C   DECLARACION DE VARIABLES

        IMPLICIT NONE

C   VARIABLES INTERNAS

        INTEGER*4 T, TI

C   VARIABLES EXTERNAS

        INTEGER*4 HY, HZ, NY, NZ, I2, S

        REAL*8 K(NY, NZ), E(NY, NZ), X1(S), Y1(S),
     $         Z(NZ), U0(NZ), K0(NZ), E0(NZ), HM, Z2(S)

C   COMIENZO DE SUBRUTINA


C        IF (I2.EQ.S) THEN
C          IF (HM.LE.Z(1)) TI=1
C           IF (HM.GE.Z(NZ)) TI=NZ-1
C           DO T=1, NZ-1, 1
C              IF ((Z(T).LE.HM).AND.(HM.LE.Z(T+1))) TI=T
C           ENDDO
C           WRITE(90,FMT='(A,I3,2(A,F9.0),5(A,F9.4))')
C     $               ' MAST    ',I2,' X = ',X1(I2),' Y = ',Y1(I2),
C     $               ' K = ',K(HY, HZ),' E = ',E(HY, HZ),
C     $    ' U0 = ',U0(TI)+(HM-Z(TI))*(U0(TI+1)-U0(TI))/(Z(TI+1)-Z(TI)),
C     $    ' K0 = ',K0(TI)+(HM-Z(TI))*(K0(TI+1)-K0(TI))/(Z(TI+1)-Z(TI)),
C     $    ' E0 = ',E0(TI)+(HM-Z(TI))*(E0(TI+1)-E0(TI))/(Z(TI+1)-Z(TI))
C           ELSE
              IF (Z2(I2).LE.Z(1)) TI=1
              IF (Z2(I2).GE.Z(NZ)) TI=NZ-1
              DO T=1, NZ-1, 1
                 IF ((Z(T).LE.Z2(I2)).AND.(Z2(I2).LE.Z(T+1))) TI=T
              ENDDO
              WRITE(90,FMT='(A,I3,2(A,F9.0),5(A,F9.4))')
     $               ' TURBINE ',I2,' X = ',X1(I2),' Y = ',Y1(I2),
     $               ' K = ',K(HY, HZ),' E = ',E(HY, HZ),
     $ ' U0 = ',U0(TI)+(Z2(I2)-Z(TI))*(U0(TI+1)-U0(TI))/(Z(TI+1)-Z(TI)),
     $ ' K0 = ',K0(TI)+(Z2(I2)-Z(TI))*(K0(TI+1)-K0(TI))/(Z(TI+1)-Z(TI)),
     $ ' E0 = ',E0(TI)+(Z2(I2)-Z(TI))*(E0(TI+1)-E0(TI))/(Z(TI+1)-Z(TI))
C        ENDIF

        RETURN
        END

C ==========================================================================
