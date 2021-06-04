C****************************************************************************
C
C		 	EPS09.f
C
C An interface for the scale dependent nuclear modifications
C 		R_f^A(x,Q) = f_A(x,Q)/f_p(x,Q) 
C where f_A is the distribution of the parton flavour f for a PROTON in a
C nucleus A, and f_p is the corresponding parton distribution in the 
C free proton.
C  
C When using this interface, please refer to:
C  
C K.J. Eskola, H. Paukkunen and C.A. Salgado,
C "EPS09 - a New Generation of NLO and LO Nuclear Parton Distribution Functions,"
C Published as JHEP04(2009) 065.
C Eprint: arXiv:0902.4154 [hep-ph].
C
C Questions & comments to:
C   hannu.paukkunen@phys.jyu.fi
C   kari.eskola@phys.jyu.fi
C   carlos.salgado@usc.es
C 
C ***************************************************************************
C Instructions:
C
C For given input values of
C
C     order: 1=LO, 2=NLO   ; integer
C     pset : 1...31        ; integer
C            1     = central fit
C            2,3   = error sets S{+1}, S{-1}
C            4,5   = error sets S{+2}, S{-2}
C            ...   ...
C            30,31 = error sets {S+15}, {S-15}
C     A    : atomic number ; integer
C     x    : Bjorken-x     ; double precision
C     Q    : scale in GeV  ; double precision
C
C the command 
C
C   Call EPS09(order, pset, A, x, Q, ruv, rdv, ru, rd, rs, rc, rb, rg)
C
C returns the bound proton nuclear corrections R_f^A(x,Q)
C (in double precision) for
C	
C	ruv = up valence
C	rdv = down valence
C	ru  = up sea
C	rd  = down sea
C	rs  = strange
C	rc  = charm
C	rb  = bottom
C	rg  = gluons
C
C The nuclear corrections for bound neutrons can be obtained
C by the isospin symmetry, e.g. the total up quark distribution
C per nucleon in a nucleus A with Z protons is
C
C  u_A(x,Q) =    Z/A * [ruv*uV_p(x,Q) + ru*uSea_p(x,Q)] +
C            (A-Z)/A * [rdv*dV_p(x,Q) + rd*dSea_p(x,Q)]
C
C Note that the parametrization should only be applied at the
C kinematical domain
C
C             1e-6 <= x <= 1
C              1.3 <= Q <= 1000 GeV.
C
C No warning message is displayed if these limits are
C exceeded, and outside these boundaries the modifications
C are frozen to the boundary values, i.e
C
C   for Q > 1000, the modifications at Q=1000 are returned,
C   for Q < 1.3,  the modifications at Q=1.3 are returned,
C   for x > 1,    the modifications at x=1 are returned
C   for x < 1e-6, the modifications at x=1e-6 are returned,
C
C The data used by the program for required order
C and atomic number A, are stored in separate files
C
C   LO : EPS09LOR_A
C   NLO: EPS09NLOR_A
C
C which must be located in the working directory.
C
C The error bands for absolute cross-sections and for
C their nuclear ratios should be computed as explained
C in Secs. 2.5 and 4 of arXiv:0902.4154 [hep-ph]. For
C the absolute cross sections, both the errors in the
C free-proton PDFs f_p(x,Q) and the errors in
C the modifications R_f^A(x,Q) should be accounted for.
C For the nuclear ratios, it is sufficient to account only
C for the errors in the modifications R_f^A(x,Q).
C
C *********************************************************
C *********************************************************
      program main
      implicit none
      integer order, pset, A
      integer k
      Double precision :: ruv, rdv, ru, rd, rs, rc, rb, rg
      Double precision :: Q, x

      Double precision :: pruv, prdv, pru, prd, prs, prc, 
     &                    prb, prg
      Double precision :: mruv, mrdv, mru, mrd, mrs, mrc, 
     &                    mrb, mrg
      Double precision :: p2ruv, p2rdv, p2ru, p2rd, p2rs, p2rc, 
     &                    p2rb, p2rg
      Double precision :: m2ruv, m2rdv, m2ru, m2rd, m2rs, m2rc, 
     &                    m2rb, m2rg
      Double precision :: p3ruv, p3rdv, p3ru, p3rd, p3rs, p3rc, 
     &                    p3rb, p3rg
      Double precision :: m3ruv, m3rdv, m3ru, m3rd, m3rs, m3rc, 
     &                    m3rb, m3rg

      Double precision :: p4ruv, p4rdv, p4ru, p4rd, p4rs, p4rc, 
     &                    p4rb, p4rg
      Double precision :: m4ruv, m4rdv, m4ru, m4rd, m4rs, m4rc, 
     &                    m4rb, m4rg
      Double precision :: p5ruv, p5rdv, p5ru, p5rd, p5rs, p5rc, 
     &                    p5rb, p5rg
      Double precision :: m5ruv, m5rdv, m5ru, m5rd, m5rs, m5rc, 
     &                    m5rb, m5rg
      Double precision :: p6ruv, p6rdv, p6ru, p6rd, p6rs, p6rc, 
     &                    p6rb, p6rg
      Double precision :: m6ruv, m6rdv, m6ru, m6rd, m6rs, m6rc, 
     &                    m6rb, m6rg
      Double precision :: p7ruv, p7rdv, p7ru, p7rd, p7rs, p7rc, 
     &                    p7rb, p7rg
      Double precision :: m7ruv, m7rdv, m7ru, m7rd, m7rs, m7rc, 
     &                    m7rb, m7rg

      Double precision :: p8ruv, p8rdv, p8ru, p8rd, p8rs, p8rc, 
     &                    p8rb, p8rg
      Double precision :: m8ruv, m8rdv, m8ru, m8rd, m8rs, m8rc, 
     &                    m8rb, m8rg
      Double precision :: p9ruv, p9rdv, p9ru, p9rd, p9rs, p9rc, 
     &                    p9rb, p9rg
      Double precision :: m9ruv, m9rdv, m9ru, m9rd, m9rs, m9rc, 
     &                    m9rb, m9rg
      Double precision :: p10ruv, p10rdv, p10ru, p10rd, p10rs, p10rc, 
     &                    p10rb, p10rg
      Double precision :: m10ruv, m10rdv, m10ru, m10rd, m10rs, m10rc, 
     &                    m10rb, m10rg
      Double precision :: p11ruv, p11rdv, p11ru, p11rd, p11rs, p11rc, 
     &                    p11rb, p11rg
      Double precision :: m11ruv, m11rdv, m11ru, m11rd, m11rs, m11rc, 
     &                    m11rb, m11rg

      Double precision :: p12ruv, p12rdv, p12ru, p12rd, p12rs, p12rc, 
     &                    p12rb, p12rg
      Double precision :: m12ruv, m12rdv, m12ru, m12rd, m12rs, m12rc, 
     &                    m12rb, m12rg
      Double precision :: p13ruv, p13rdv, p13ru, p13rd, p13rs, p13rc, 
     &                    p13rb, p13rg
      Double precision :: m13ruv, m13rdv, m13ru, m13rd, m13rs, m13rc, 
     &                    m13rb, m13rg
      Double precision :: p14ruv, p14rdv, p14ru, p14rd, p14rs, p14rc, 
     &                    p14rb, p14rg
      Double precision :: m14ruv, m14rdv, m14ru, m14rd, m14rs, m14rc, 
     &                    m14rb, m14rg
      Double precision :: p15ruv, p15rdv, p15ru, p15rd, p15rs, p15rc, 
     &                    p15rb, p15rg
      Double precision :: m15ruv, m15rdv, m15ru, m15rd, m15rs, m15rc, 
     &                    m15rb, m15rg
      Double precision :: p16ruv, p16rdv, p16ru, p16rd, p16rs, p16rc, 
     &                    p16rb, p16rg
      Double precision :: m16ruv, m16rdv, m16ru, m16rd, m16rs, m16rc, 
     &                    m16rb, m16rg
     
      Double precision :: druv, drdv, dru, drd, drs, drc, drb, drg
      Double precision :: pT, rootS 

c      Write(*,*) 'Executing main'
      order = 1
      pset = 1
      Q = 10.0d0
      x = 0.1d0
      A = 208
      rootS = 5020.0d0

      do k = 1, 30, 1 
           pT = 100.d0 + 25.0d0*(k-1)
           x = 2*pT/rootS
           Call EPS09(order, pset, A, x, Q, ruv, rdv, ru, rd, rs, rc, 
     &          rb, rg)
           Call EPS09(order, 2, A, x, Q, 
     &          pruv, prdv, pru, prd, prs, prc, 
     &          prb, prg)
           Call EPS09(order, 3, A, x, Q, 
     &          mruv, mrdv, mru, mrd, mrs, mrc, 
     &          mrb, mrg)
           Call EPS09(order, 4, A, x, Q, 
     &          p2ruv, p2rdv, p2ru, p2rd, p2rs, p2rc, 
     &          p2rb, p2rg)
           Call EPS09(order, 5, A, x, Q, 
     &          m2ruv, m2rdv, m2ru, m2rd, m2rs, m2rc, 
     &          m2rb, m2rg)
           Call EPS09(order, 6, A, x, Q, 
     &          p3ruv, p3rdv, p3ru, p3rd, p3rs, p3rc, 
     &          p3rb, p3rg)
           Call EPS09(order, 7, A, x, Q, 
     &          m3ruv, m3rdv, m3ru, m3rd, m3rs, m3rc, 
     &          m3rb, m3rg)
           Call EPS09(order, 8, A, x, Q, 
     &          p4ruv, p4rdv, p4ru, p4rd, p4rs, p4rc, 
     &          p4rb, p4rg)
           Call EPS09(order, 9, A, x, Q, 
     &          m4ruv, m4rdv, m4ru, m4rd, m4rs, m4rc, 
     &          m4rb, m4rg)

           Call EPS09(order, 10, A, x, Q, 
     &          p5ruv, p5rdv, p5ru, p5rd, p5rs, p5rc, 
     &          p5rb, p5rg)
           Call EPS09(order, 11, A, x, Q, 
     &          m5ruv, m5rdv, m5ru, m5rd, m5rs, m5rc, 
     &          m5rb, m5rg)
           Call EPS09(order, 12, A, x, Q, 
     &          p6ruv, p6rdv, p6ru, p6rd, p6rs, p6rc, 
     &          p6rb, p6rg)
           Call EPS09(order, 13, A, x, Q, 
     &          m6ruv, m6rdv, m6ru, m6rd, m6rs, m6rc, 
     &          m6rb, m6rg)
           Call EPS09(order, 14, A, x, Q, 
     &          p7ruv, p7rdv, p7ru, p7rd, p7rs, p7rc, 
     &          p7rb, p7rg)
           Call EPS09(order, 15, A, x, Q, 
     &          m7ruv, m7rdv, m7ru, m7rd, m7rs, m7rc, 
     &          m7rb, m7rg)
           Call EPS09(order, 16, A, x, Q, 
     &          p8ruv, p8rdv, p8ru, p8rd, p8rs, p8rc, 
     &          p8rb, p8rg)
           Call EPS09(order, 17, A, x, Q, 
     &          m8ruv, m8rdv, m8ru, m8rd, m8rs, m8rc, 
     &          m8rb, m8rg)

           Call EPS09(order, 18, A, x, Q, 
     &          p9ruv, p9rdv, p9ru, p9rd, p9rs, p9rc, 
     &          p9rb, p9rg)
           Call EPS09(order, 19, A, x, Q, 
     &          m9ruv, m9rdv, m9ru, m9rd, m9rs, m9rc, 
     &          m9rb, m9rg)
           Call EPS09(order, 20, A, x, Q, 
     &          p10ruv, p10rdv, p10ru, p10rd, p10rs, p10rc, 
     &          p10rb, p10rg)
           Call EPS09(order, 21, A, x, Q, 
     &          m10ruv, m10rdv, m10ru, m10rd, m10rs, m10rc, 
     &          m10rb, m10rg)
           Call EPS09(order, 22, A, x, Q, 
     &          p11ruv, p11rdv, p11ru, p11rd, p11rs, p11rc, 
     &          p11rb, p11rg)
           Call EPS09(order, 23, A, x, Q, 
     &          m11ruv, m11rdv, m11ru, m11rd, m11rs, m11rc, 
     &          m11rb, m11rg)
           Call EPS09(order, 24, A, x, Q, 
     &          p12ruv, p12rdv, p12ru, p12rd, p12rs, p12rc, 
     &          p12rb, p12rg)
           Call EPS09(order, 25, A, x, Q, 
     &          m12ruv, m12rdv, m12ru, m12rd, m12rs, m12rc, 
     &          m12rb, m12rg)

           Call EPS09(order, 26, A, x, Q, 
     &          p13ruv, p13rdv, p13ru, p13rd, p13rs, p13rc, 
     &          p13rb, p13rg)
           Call EPS09(order, 27, A, x, Q, 
     &          m13ruv, m13rdv, m13ru, m13rd, m13rs, m13rc, 
     &          m13rb, m13rg)
           Call EPS09(order, 28, A, x, Q, 
     &          p14ruv, p14rdv, p14ru, p14rd, p14rs, p14rc, 
     &          p14rb, p14rg)
           Call EPS09(order, 29, A, x, Q, 
     &          m14ruv, m14rdv, m14ru, m14rd, m14rs, m14rc, 
     &          m14rb, m14rg)
           Call EPS09(order, 30, A, x, Q, 
     &          p15ruv, p15rdv, p15ru, p15rd, p15rs, p15rc, 
     &          p15rb, p15rg)
           Call EPS09(order, 31, A, x, Q, 
     &          m15ruv, m15rdv, m15ru, m15rd, m15rs, m15rc, 
     &          m15rb, m15rg)
           Call EPS09(order, 24, A, x, Q, 
     &          p16ruv, p16rdv, p16ru, p16rd, p16rs, p16rc, 
     &          p16rb, p16rg)
           Call EPS09(order, 25, A, x, Q, 
     &          m16ruv, m16rdv, m16ru, m16rd, m16rs, m16rc, 
     &          m16rb, m16rg)


c           drg = sqrt((prg-mrg)*(prg-mrg)/4.0)
           drg = sqrt( (prg-mrg)*(prg-mrg)/4.0+
     &                 (p2rg-m2rg)*(p2rg-m2rg)/4.0+
     &                 (p3rg-m3rg)*(p3rg-m3rg)/4.0+
     &                 (p4rg-m4rg)*(p4rg-m4rg)/4.0+
     &                 (p5rg-m5rg)*(p5rg-m5rg)/4.0+
     &                 (p6rg-m6rg)*(p6rg-m6rg)/4.0+
     &                 (p7rg-m7rg)*(p7rg-m7rg)/4.0+
     &                 (p8rg-m8rg)*(p8rg-m8rg)/4.0+
     &                 (p9rg-m9rg)*(p9rg-m9rg)/4.0+
     &                 (p10rg-m10rg)*(p10rg-m10rg)/4.0+
     &                 (p11rg-m11rg)*(p11rg-m11rg)/4.0+
     &                 (p12rg-m12rg)*(p12rg-m12rg)/4.0+
     &                 (p13rg-m13rg)*(p13rg-m13rg)/4.0+
     &                 (p14rg-m14rg)*(p14rg-m14rg)/4.0+
     &                 (p15rg-m15rg)*(p15rg-m15rg)/4.0+
     &                 (p16rg-m16rg)*(p16rg-m16rg)/4.0   )
c           Write(*,*) pT, x, rg, prg, mrg
            Write(*,*) pT, x, rg, drg
      end do
      end program

      Subroutine EPS09(order, pset, AAA, xxx, QQQ,
     &                   ruv, rdv, ru, rd, rs, rc, rb, rg)

      Implicit none
      Double precision :: ruv, rdv, ru, rd, rs, rc, rb, rg, QQQ, xxx
      Double precision :: LSTEP, x, Q, Q2, allvalues(1:31,1:8,0:50,0:50)
      Double precision :: x_i=0.000001, arg(4), fu(4), res, fg(3)
      Double precision :: result(9), dummy
      Double precision :: realQ, Q2min=1.69, Q2max=1000000.0, Qsteps=50.0
      Double precision :: n_x, zero=0.0

      Character (Len=50) filenimi

      Integer :: xlinsteps=25, xlogsteps=25, startline, lineno
      Integer :: k, p, t, Qpoint, xpoint, pset, iovar
      Integer :: setnumber,j, A, openchannel, order, AAA
      Integer :: psetlast = -10, Alast = -10, orderlast = -10

      save Alast
      save psetlast
      save orderlast
      save allvalues

C *********************************************
C Stop if the set specifications are wrong ones
C *********************************************

      If (order .NE. 1 .and. order .NE. 2) then
      Write(*,*) 'Wrong order!'
      Write(*,*) 'LO : order = 1'
      Write(*,*) 'NLO: order = 2'
      Stop
      End If

      If (pset .LT. 1 .or. pset .GT. 31) then
      Write(*,*) 'Wrong set!'
      Write(*,*) 'Central set: pset = 1'
      Write(*,*) 'Error sets : pset = 2...31'
      Stop
      End If

C ********************************
C Make sure not to change any
C specifications given by the user
C ********************************

      A  = AAA
      x  = xxx
      Q  = QQQ
      Q2 = Q*Q 

C *******************************
C Freeze x if it's < 10E-6 or > 1
C *******************************

      If (x .LT. x_i) Then
      x = x_i
      End If
      If (x .GT. 1) Then
      x = 1.0
      End If

C ************************************
C Freeze Q^2 if it's < 1.69 or > 10E+6
C ************************************

      If (Q2 .LT. Q2min) Then
      Q2 = Q2min
      End If
      If (Q2 .GT. Q2max) Then
      Q2 = Q2max
      End If

C If the set specifications have been changed, read the tables again

      If (A .NE. Alast .or.
     &   order .NE. orderlast) Then

C      Write(*,*) 'Set changed!'

C Read the table

      If (order .EQ. 1) then

        If (A < 10) Then
        Write(filenimi,'("EPS09LOR_", I1)'), A
        Else If (A < 100) Then
        Write(filenimi,'("EPS09LOR_", I2)'), A
        Else If (A < 1000) Then
        Write(filenimi,'("EPS09LOR_", I3)'), A
        End If

      Else

        If (A < 10) Then
        Write(filenimi,'("EPS09NLOR_", I1)'), A
        Else If (A < 100) Then
        Write(filenimi,'("EPS09NLOR_", I2)'), A
        Else If (A < 1000) Then
        Write(filenimi,'("EPS09NLOR_", I3)'), A
        End If

      End If

      Call NextUnit(openchannel)

      OPEN (openchannel, file = filenimi, status='OLD', IOSTAT=iovar)

      If (iovar .NE. 0) Then
      Write(*,*) 'Missing file: ',filenimi
      stop
      End If

      Do setnumber = 1, 31

      Do k = 0,50

      Read(openchannel,*) dummy

      Do t = 0,50

      Read(openchannel,*) (allvalues(setnumber,p,k,t), p=1,8)

      End Do
      End Do

      End Do

      Close(openchannel)

      psetlast  = pset
      Alast     = A
      orderlast = order

      End If

C Find out the position in the loglog Q^2-grid

      realQ  = Qsteps * (log(log(Q2)/log(Q2min)))/
     &                  (log(log(Q2max)/log(Q2min)))
      Qpoint = Aint(realQ)

      If (Qpoint .LE. 0) Then
         Qpoint = 1
      End If
      If (Qpoint .GE. Anint(Qsteps)-1) Then
         Qpoint = Anint(Qsteps)-1
      End If

      LSTEP = (1.0/(xlogsteps)) * LOG(0.1/x_i)

C *********************
C Interpolate the grids 
C *********************

      Do t=1,8

C Find the position in the x-grid

      If (x .LE. 0.1) then
         n_x  = ((1.0/LSTEP) * Log(x/x_i))
       xpoint = Aint(n_x)
      Else
       n_x    = ((x-0.1)*xlinsteps/(1.0-0.1) + xlogsteps)
       xpoint = Aint(n_x)
      End If

      If (xpoint .LE. 0) Then
        xpoint = 1
      End If

      If (t .EQ. 1 .or. t .EQ. 2) Then
         If (xpoint .GE. (xlinsteps+xlogsteps)-4) Then
         xpoint = (xlinsteps+xlogsteps)-4
         End If
      End If

      If (t .EQ. 3 .or. t .EQ. 4 .or. t .EQ. 5 .or.
     &    t .EQ. 6 .or. t .EQ. 7) Then
        If (xpoint .GE. (xlinsteps+xlogsteps)-7) Then
        xpoint = (xlinsteps+xlogsteps)-7
        End If
      End If
	
      If (t .EQ. 8) Then
        If (xpoint .GE. (xlinsteps+xlogsteps)-4) Then
        xpoint = (xlinsteps+xlogsteps)-4
        End If
      End If

      Do k = 1, 4
      If (xpoint-2+k .LT. xlogsteps) Then
      arg(k) = (x_i) * exp(LSTEP * (xpoint-2+k))
      Else
      arg(k) = 0.1 + (xpoint-2+k-xlogsteps) * (1-0.1)/xlinsteps
      End If
      End Do
 
      Do j=1,3

      fu(1) = allvalues(pset,t,Qpoint-2+j,xpoint-1)
      fu(2) = allvalues(pset,t,Qpoint-2+j,xpoint)
      fu(3) = allvalues(pset,t,Qpoint-2+j,xpoint+1)
      fu(4) = allvalues(pset,t,Qpoint-2+j,xpoint+2)
      Call luovi(fu,arg,4,x,res)
      fg(j) = res

C *****************************************
C *****************************************

      End Do

      arg(1) = Qpoint-1
      arg(2) = Qpoint
      arg(3) = Qpoint+1

      Call luovi(fg,arg,3,realQ,res)
  
      result(t) = res

      End Do

      ruv = max(result(1),zero)
      rdv = max(result(2),zero)
      ru  = max(result(3),zero)
      rd  = max(result(4),zero)
      rs  = max(result(5),zero)
      rc  = max(result(6),zero)
      rb  = max(result(7),zero)
      rg  = max(result(8),zero)

200   Continue

      End Subroutine EPS09

! ********************************
! Modified version of Cern Library
! interpolation routine E100
! ********************************

      SUBROUTINE luovi(F,ARG,MMM,Z,SUM)

      Implicit none
      INTEGER :: MMM
      Double precision  :: F(MMM), ARG(MMM), COF(MMM), SUM, Z
      INTEGER :: M, MM, I, J, JNDEX, INDEX

      MM = MIN(MMM, 20)
      M = MM - 1
      DO 1780 I= 1, MM
      COF(I) = F(I)
 1780 Continue
      DO 1800 I= 1, M
      DO 1790 J= I, M
         JNDEX = MM - J
         INDEX = JNDEX + I
         COF(INDEX) = (COF(INDEX)-COF(INDEX-1))/(ARG(INDEX)-ARG(JNDEX))
 1790 CONTINUE
 1800 CONTINUE
      SUM = COF(MM)
      DO 1810 I= 1, M
         INDEX = MM - I
         SUM = (Z-ARG(INDEX))*SUM + COF(INDEX)
 1810 CONTINUE

      End SUBROUTINE luovi

! **********************
! Find the open I/O unit
! **********************

      Subroutine NextUnit(firstopen)

      Logical EX
      Integer firstopen, N

      Do 10 N = 10, 300
C        write(*,*) N
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            firstopen = N
            Goto 20 
        Endif
10    Continue
      Stop ' There is no available I/O unit. '
20    Continue
      End Subroutine NextUnit
