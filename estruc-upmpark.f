

C     *******************************************************
C-4-  PROGRAM UPMWAKE TO CALCULATE WIND-TURBINE WAKES
C     *******************************************************
C-7-    INPUT DATA
C       ==========
C    ****************************************
C-81-       COMIENZO DEL PROGRAMA
C    ****************************************
C-85- DECLARACION DE VARIABLES
C-89- PARAMETROS
C-108- DICCIONARIO DE PARAMETROS UTILES PARA EL USUARIO
C-140- VARIABLES LOGICAS
C-144- VARIABLES REALES
C-191- VARIABLES ENTERAS
C-205- VARIABLES TIPO CARACTER
    write(178,*)randomreal(I)
C-239 ***** INTRODUCCION DEL NOMBRE DEL FICHERO DE DATOS DEL PARQUE *****
1       WRITE (*, *) 'PARK DATA FILE? '
 WRITE(*,*) DATA$
    WRITE (*,*) 'REMEMBER ".dat"'
C-258- ****** LECTURA FICHERO INICIALIZACION DEL PARQUE *****
       WRITE(*,fmt='(2a)')'Programa no preparado para ',
    WRITE (*,*) 'NO SE SELECCIONA UN TIPO CORRECTO DE MEANDERING'
       WRITE (*,*) 'NO SE SELECCIONA UN TIPO CORRECTO DE SPECTRO'
          WRITE(*,FMT='(3A,F8.3,A,F8.3,A)')'EL L�MITE SUPERIOR',
          WRITE(*,FMT='(3A,F8.3,A)')'POR TANTO, EL L�MITE ',
       WRITE (*,*) 'TIPO DE PONDERACI�N INVALIDO'
       WRITE (*,*) 'ENTRADA DE KMWM NO V�LIDA'
             WRITE(*,FMT='(3a)')" ESTUDIO DEL CASO 1 DEL KMWM:
             WRITE(*,FMT='(3a)')" ESTUDIO DEL CASO 2 DEL KMWM:
             WRITE(*,FMT='(5a)')" ESTUDIO DEL CASO 3 DEL KMWM:

    WRITE (93,'(A)')DATA$
       WRITE(95,fmt='(2a)')' VINICIAL SPECTRUMX SPECTRUMY ',
C  ****************************
C-432-  *** BUCLE DE DIRECCIONES ***
C  ****************************
C  *********************************************************************************************
C-469-  *** MODIFICACION DE LA VELOCIDAD Y DE LA DIRECCI�N PARA EL CASO DE ESTUDIO DEL MEANDERING ***
C  *********************************************************************************************
    WRITE(95,fmt='(9f12.5)')VHMEMORY(A1),SPECTRUMX,SPECTRUMY,
C-517- ***** APERTURA DE FICHEROS DE SALIDA ******
    WRITE(AUXILIAR9,FMT='(a1,I6.6,a2)')'t',ctiempo,'cs'
C-578- ***** DICCIONARIO DE NUMEROS ASOCIADOS A CADA FICHERO *****
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
!              WRITE (94, *) (SECMAT(SEC, J), J = 1 , 9)
!           WRITE (94, *) (SECMAT(1, J), J = 1 , 9)
!        WRITE (94, '(I3)') S
!           WRITE (94, '(F9.1)') X1(I)
!           WRITE (94, '(F9.1)') Y1(I)
!           WRITE (94, '(F8.3)') Z2(I)
!           WRITE (94, '(A16)') TURB$(I)
!        WRITE (94, '(F8.3)') HM
!        WRITE (94, '(I3)') ND
!           WRITE (94, '(F8.3)') DI(I)
!           WRITE (94, '(F8.3)') VH(I)
!        WRITE (94,'(I3)') NTS
!           WRITE (94, '(F8.3)') XTS(I)
!           WRITE (94, '(F8.3)') YTS(I)
!           WRITE (94, '(F8.3)') ZTS(I)
!        WRITE (94,'(F8.3)') TDWA
!        WRITE (94, '(I3)') CASO
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
!              WRITE(94,'(I3)') NSKMWM
!                 WRITE(94,'(I3)') SECKMWM(I)
       WRITE(*,fmt='(4a)')' Error. El fichero del tipo de ',
C-755- ***** OROGRAFIA: LECTURA DEL FICHERO DE OROGRAFIA *****
C-766- ***** LECTURA DE DIAMETROS *****
         write(*,fmt='(a,i2.2,a)')' El di�metro de la turbina ',
         write(*,fmt='(a,i2.2,a)')' Error1: La unidad ',
C-785- ****** CALCULO DE LA VELOCIDAD INCIDENTE CORREGIDA ********
C-793- ****** ALTURA DEL MASTIL METEOROLOGICO (SI EXISTE) *****
C-799- ****** ENCABEZAMIENTO PARA CADA CASO DE ANALISIS ******
C-849- ****** CAMBIO DE EJES ******
C *********************************************
C-857- *******  CREACION DE LA MALLA INICIAL *******
C *********************************************
C-861- ***** CALCULO DE INCREMENTOS DE LA MALLA DY, DZ *****
C-866- ***** DETERMINACION DE LOS LIMITES DE LA MALLA *****
       WRITE(93,'(A)')'Nombre de ficheros de casos calculados'
    WRITE(93,FMT='(A)')DATA$(:R) // AUXILIAR9
     WRITE (*, 10001) NY, NZ
C     ****************************
C-901-     ******* FLUJO BASICO *******
C     ****************************
C-906-     CALCULO DEL FLUJO BASICO
C-912-     IMPRESION FICHERO DEL FLUJO BASICO
C     *************************************
C-919-     ***** CALCULOS EN CADA SECCION ******
C     *************************************
C-923-     INICIALIZACION DE LOS NUDOS DE LA MALLA
C-929-     COMIENZO DEL BUCLE DE SECCIONES
 WRITE (*, 10002) NP
42      WRITE (*, 10003) N
C ***********************************************************
C-957- ********* ANALISIS DEL ESPECTRO DE LA TURBULENCIA *********
C ***********************************************************
C-969- ********* ANALISIS DE INTERACCION ***********
C *********************************************
C **************************************
C-981- ********REDEFINICION DE MALLA*********
C **************************************
   WRITE (*, 10004) I
C-992-       DETERMINACION DE LA POSICION DEL ROTOR DE INTERACCION EN LA MALLA
C-997-       DETERMINACION DE LA VELOCIDAD EN EL CENTRO DEL ROTOR
C-1001-       DETERMINACION DE LA K EN EL CENTRO DEL ROTOR
C-1014-       LECTURA DE LAS CARACTERISTICAS DE LA TURBINA INTERCEPTADA
C-1019-       CAMBIO DE LAS CONDICIONES DE CONTORNO EN LA MALLA
C-1033-       EFECTO DE LA OROGRAFIA EN LA POTENCIA
C ************************************
C-1046- ********CALCULO DE POTENCIA*********
C ************************************
C **************************************
C-1054- **** ALMACENAMIENTO DE RESULTADOS ****
C **************************************
C-1063-     AVERIGUAMOS SI LA SECCION ACTUAL ESTA ESPECIFICADA EN EL FICHERO
C-1064-     DE DEFINICION DEL PARQUE PARA ALMACENAR DATOS
C ***************************
C-1123- *****�ULTIMA SECCION?******
C ***************************
C **********************
C-1135- **** RENDIMIENTOS ****
C **********************
C ********************************
C-1192- ****** �ULTIMO INSTANTE?  ******
C ********************************
C ********************************
C-1213- ****** �ULTIMA DIRECCION? ******
C ********************************
c      ===========================================================
c-1222-                        FICHERO GRAFICO
c      ===========================================================
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',gru$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csu.header'
       write(n,*) 'file = ',gru$
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grv$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csv.header'
       write(n,*) 'file = ',grv$
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grw$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csw.header'
       write(n,*) 'file = ',grw$
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grp$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csp.header'
       write(n,*) 'file = ',grp$
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grt$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'cst.header'
       write(n,*) 'file = ',grt$
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grk$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csk.header'
       write(n,*) 'file = ',grk$
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',gre$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'cse.header'
       write(n,*) 'file = ',gre$
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',gvi$
!        write(fichero18,fmt='(a1,i6.6,a11)')'t',ctiempo,'csvi.header'
       write(n,*) 'file = ',gvi$
!        write(n,*) 'file = /u/emilio/EOLICA/GRAFICA/',grm$
!        write(fichero17,fmt='(a1,i6.6,a10)')'t',ctiempo,'csm.header'
       write(n,*) 'file = ',grm$
    write(n,*) 'grid = ',np,' x ',nz,' x ',ny-2
    write(n,*) 'field = scalar_one'
    write(n,*) 'type = float'
    write(n,*) 'format = ascii'
    write(n,*) 'structure = scalar'
    write(n,*) 'interleaving = record'
    write(n,*) 'header = lines 12'
    write(n,2222)  dx,dz,dy
C ****************************************************
C-1313- ****** PONDERACI�N DE LOS CASOS DE MEANDERING ******
C ****************************************************
C ********************************
C-1325- ******* FIN DEL PROGRAMA *******
C ********************************
C **************************
C-1344- ******* SUBRUTINAS *******
C **************************
C --------------------------------------------------------------------------
C-1352-     Subrutina que determina la altura a la que se mide la velocidad in-
C  cidente en el parque en funcion de que exista o no mastil meteorologico.
C --------------------------------------------------------------------------
C  ----------------------------------------------------------------------
C-1398-     Subrutina que calcula las magnitudes que describen el flujo basico
C  del lugar en cuestion, para cada caso de analisis.
C  ----------------------------------------------------------------------
C  ----------------------------------------------------------------------
C-1517-     Subrutina que resuelve las ecuaciones diferenciales del modelo, y
C  calcula las magnitudes U, V, W, T, P, K, E, VI.
C  ----------------------------------------------------------------------
c original       WRITE (*,10000) IT, Er
       WRITE (*,10001) IT, Er, maxpp, N, NP-N, A1, ND, DATA$
        WRITE (*,10002) IT, Er, maxpp, N, NP-N, A1, ND, DATA$
C -----------------------------------------------------------
C-2156-  CALCULO DE ANALISIS DEL Kinematic Model For Wake Meandering
C -----------------------------------------------------------
                      WRITE(*,*)'NO SE EMPLEA ESTE C�LCULO'
C ==========================================================================
C  -----------------------------------------------------------
C-2216-  Subrutina para el cambio de los ejes arbitrarios iniciales
C  a otros ejes (X2,Y2,Z2), en los que el eje X2 tiene la direc-
C  cion del viento
C  -----------------------------------------------------------
C  ------------------------------------------------------------------
C-2267-    Subrutina que calcula la velocidad incidente en el rotor para
C  la estimacion de la potencia que produce la aeroturbina.
C  ------------------------------------------------------------------
C  -----------------------------------------------------------------------
C-2360-       Subrutina que corrige la velocidad medida en el mastil segun la
C   influencia de la orografia en el emplazamiento del mismo; a esta
C   velocidad corregida estan referidos el resto de los porcentajes
C   de influencia de la orografia en los emplazamientos de las turbinas.
C  -----------------------------------------------------------------------
              WRITE (*,10000) NM(A1)
3            WRITE (*,10001) AB1
      WRITE (*,*) 'NUMBER OF THE MAST THAT MEASURES WIND ',
         WRITE (*,10003) NT
C  --------------------------------------------------------------------
C-2461-     Subrutina que a�ade al fichero de potencias y rendimientos el
C  encabezamiento que separa cada caso de analisis. Tambien lo indica
C  por pantalla.
C  --------------------------------------------------------------------
     WRITE (FIC,10000)
     WRITE (FIC,10001)
     WRITE (FIC,10002)
     WRITE (FIC,10003)
     WRITE (FIC,10004)
     WRITE (FIC,10005)
     WRITE (FIC,10006)
     WRITE (FIC,10007)
     WRITE (FIC,10020)
     WRITE (FIC, 9998)
     WRITE (FIC, 9999)
     WRITE (FIC, 10008) DATA$
 WRITE (FIC, 10009)
    WRITE (FIC, 10010) NM1
    WRITE (FIC, 10011) ALT
    WRITE (FIC, 10012) VH(A1)
    WRITE (FIC, 10013) AB1
    WRITE (FIC, 10050) TIEMPO
       WRITE(FIC, 23456)
    WRITE (FIC, 10014) ALT
    WRITE (FIC, 10015) VH(A1)
    WRITE (FIC, 10016) AB1
    WRITE (FIC, 10051) TIEMPO
       WRITE(FIC, 12345)
 WRITE (FIC, 10017)
    WRITE (FIC, 10018)
    WRITE (FIC, 10019)
C ----------------------------------------------------------------
C-2605-     Funcion que calcula el espectro de la turbulencia segun la
C energia cinetica turbulenta y la frecuencia adimensionalizada.
C ----------------------------------------------------------------
C  ----------------------------------------------------------------
C-2634-       Subrutina que calcula el espectro de la turbulencia en los
C   puntos especificados en el fichero de definicion del parque.
C  ----------------------------------------------------------------
       WRITE (15, 10002) XTS(I), YTS(I), ZTS(I),(TT(J), J= 1, 5)
C-2695-    Subrutina que lee el fichero de orografia, formado a partir de
C los datos del codigo WASP. Procesa los datos correspondientes.
C -----------------------------------------------------------------
C-2766-  ALMACENAMIENTO DE LOS RESULTADOS DE KMWM
C  -----------------------------------------------------------------------
       WRITE(100,FMT='(2A)')' COMPROBACI�N DE KINEMATIC MODEL ',
       WRITE(100,FMT='(2A,I3.3)')' CALCULO DE DELTAK/DELTAV0^2 ',
       WRITE (100,*)
       WRITE (100, 20004)' DELTAV0: ',DELTAV0(IJ)
       WRITE (100,*)
       WRITE (100,*)'DISTANCIA EN Y'
       WRITE (100, 20001)
       WRITE (100,*)
       WRITE (100,*)'DELTAVM'
          WRITE (100, 20002)
       WRITE (100,*)
       WRITE (100,*)'DELTAK'
          WRITE (100, 20003)
       WRITE (100,*)
       WRITE (100,*)'DELTAK/DELTAV0^2'
          WRITE (100, 20000)
       WRITE (100,*)
       WRITE (100,FMT='(2A)')" **********************************
C  -----------------------------------------------------------------------
C-2845-      Subrutina que almacena en el fichero '.BSF' las magnitudes de velo-
C  cidad V, energia cinetica turbulenta K, disipacion de la energia cine-
C  tica turbulenta E y viscosidad turbulenta VI, correspondientes al flujo
C  basico y calculadas en el programa principal.
C  -----------------------------------------------------------------------
 WRITE (2, 10008)
    WRITE (2, 10009) I,Z(I),U0(I),K0(I),E0(I),VI0(I)
C  ----------------------------------------------------------------
C-2886-     Subrutina que determina los incrementos DY, DZ de la malla de
C  calculo del parque.
C  ----------------------------------------------------------------
C  ----------------------------------------------------------------
C-2953-     Subrutina que determina los incrementos DY, DZ de la malla de
C  calculo del parque.
C  ----------------------------------------------------------------
C  -------------------------------------------------------------------
C-3054-      Subrutina que calcula los limites de la malla de calculo y el
C  numero de nodos en los ejes X (NP), Y (NY) y Z (NZ).
C  -------------------------------------------------------------------
C ----------------------------------------------------------------------
C-3153-      Subrutina que promedia los distintos casos de meandering.
C ----------------------------------------------------------------------
             write(*,fmt='(2a)')
             write(*,fmt='(a,f8.2,a,f8.2,a)')
             write(*,fmt='(a)')'�Quiere continuar? (1=si, 0=no)'
       write(*,*)
       write(*,*)
       write(*,*)
    write(formato10,fmt='(a1,i4,a5)')'(',my,'f7.2)'
       write(*,fmt='(a,i3.3,a,i3.3,2a)')'Procesando el fichero ',
          write(4,fmt='(a27)')a27
          write(4,fmt='(a2)')a2
          write(4,fmt='(a2)')a2
          write(4,fmt='(a45)')a45
          write(4,fmt='(a2)')a2
          write(4,fmt='(a2)')a2
          write(4,fmt='(a58)')a58
          write(4,fmt='(a27)')a27
          write(4,fmt='(a57)')a57
          write(4,fmt='(a46)')a46
          write(4,fmt='(a19)')a19
          write(4,fmt='(a58)')a58
          write(4,fmt='(a1)')' '

          write(5,fmt='(a27)')a27
          write(5,fmt='(a2)')a2
          write(5,fmt='(a2)')a2
          write(5,fmt='(a45)')a45
          write(5,fmt='(a2)')a2
          write(5,fmt='(a2)')a2
          write(5,fmt='(a58)')a58
          write(5,fmt='(a27)')a27
          write(5,fmt='(a57)')a57
          write(5,fmt='(a46)')a46
          write(5,fmt='(a19)')a19
          write(5,fmt='(a58)')a58
          write(5,fmt='(a1)')' '
    write(*,*)'ESCRIBIENDO LOS FICHEROS DE SALIDA'
          write(4,fmt=formato10)(matrizu(i,my+1-j,kk),j=1,my,1)
          write(5,fmt=formato10)(matrizk(i,my+1-j,kk),j=1,my,1)
     write(4,fmt='(a1,//)')a1
     write(4,fmt='(a1)')a1
     write(4,fmt='(a1,//)')a1
     write(4,fmt='(a42,i1.1,a4,f5.2,/)')a42,i,a4,vt(i)

C ----------------------------------------------------------------------
C-3419-      Subrutina que estima una velocidad Uorog que resumira el efecto
C la orografia en el emplazamiento de la turbina en cuestion. Se supon-
C dra que dicho efecto se superpone linealmente al de las estelas.
C ----------------------------------------------------------------------
C-3473-     Subrutina para el calculo de la potencia de cada aeroturbina,
C  segun la velocidad cubica media calculada en el programa principal
C  y la correccion de la orografia Uorog, si es que la orografia se
C  tiene en cuenta en el analisis del parque.
C  -------------------------------------------------------------------
      WRITE(FORMATO,fmt='(a1,I2.2,a9)')'(',WARD(I),'(/),f6.2)'
      WRITE (3, 10000) I, PCU(I), I, VM(I)
      WRITE (3, 10001) I, I, VM(I)
C -------------------------------------------------------------------------
C-3570-    Subrutina para el calculo de la potencia total del parque y de los
C rendimientos de las aeroturbinas y del parque; para su determinacion
C consideraremos que la potencia maxima producida por una aeroturbina en
C el caso de analisis en cuestion corresponde a la velocidad cubica media
C en el rotor de la misma segun el flujo basico.
C -------------------------------------------------------------------------
 WRITE (3, 10000) P/1000.
          write(*,fmt='(a,i2.2,a)')' Error2: La unidad ',
       WRITE (3, 10001) I
    WRITE (3, 10003) I, RT
    WRITE (3, 10004) RP
C ----------------------------------------------------------------------
C-3724-   Subrutina de almacenamiento de resultados segun lo especificado en
C  el fichero de definicion del parque.
C ----------------------------------------------------------------------
c           WRITE (FIC, 10010) N
c           WRITE (FIC, 10011) N
c           WRITE (FIC, 10012) N
c           WRITE (FIC, 10013) N
c           WRITE (FIC, 10014) N
c           WRITE (FIC, 10015) N
c           WRITE (FIC, 10016) N
c           WRITE (FIC, 10117) N
c           WRITE (FIC, 10118) N
    WRITE (FIC, 10000) (Q(J, I), J = (NY - 1) , 2 , -1)
c        WRITE (FIC, 10018)
C ------------------------------------------------
C-3799-   Subrutina de calculo interno del metodo k-E.
C ------------------------------------------------
C ------------------------------------------------
C-3838-   Subrutina de calculo interno del metodo k-E.
C ------------------------------------------------
C  ----------------------------------------------------------------------
C-3878-     Subrutina que lee los datos del fichero de turbina correspondiente
C  y que interpola el coeficiente de empuje Ct.
C  ----------------------------------------------------------------------
    WRITE (*,10000)
1000      FORMAT (/,1X, 'Ct CURVE NOT AVAILABLE')
    WRITE (*,10001) I, VM1
    WRITE (*,*) ' Ct VALUE? '
      WRITE (*,10003)
      WRITE (*,10004) I, VM1, V1
      WRITE (*,*) ' Ct VALUE? '
      WRITE(*,'(A, F6.4)')'CT =', CT
         WRITE (*,10003)
         WRITE (*,10006) I, VM1, V1
         WRITE (*,*) ' Ct VALUE? '
         WRITE(*,'(A, F6.4)')'CT =', CT
               WRITE (91, *)'==================================='
               WRITE (91,10007) I, VM1, CT, DIAMCAMBIADO
               WRITE (91,10008) CT, CTMAX
               WRITE (91,10009) I, VM1, CT, DIAMCAMBIADO
               WRITE (91, *)'==================================='
               WRITE (91,10010) I, VM1, CT, DIAMCAMBIADO
C  -----------------------------------------------------------------------
C-4016-     Subrutina que calcula la velocidad incidente en el centro del rotor.
C  -----------------------------------------------------------------------
C  -----------------------------------------------------------------------
C-4047-     Subrutina que almacena la K y E incidente en el centro del rotor.
C  -----------------------------------------------------------------------
C           WRITE(90,FMT='(A,I3,2(A,F9.0),5(A,F9.4))')
       WRITE(90,FMT='(A,I3,2(A,F9.0),5(A,F9.4))')
