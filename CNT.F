      PROGRAM mk_cnt
C
C  Atomic coordinates for multiple unit cells of a (n,m) nanotube.
C
C  Based on routines by R. Saito published on 
C   "Physical Properties of Carbon Nanotubes" by
C  R. Saito, G. Dresselhaus and M. S. Dresselhaus
C  Imperial College Press 1995
C
C
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (nk=20000,aa=1.42)
      DIMENSION x(nk),y(nk),z(nk)
      DIMENSION iic(nk),ic(nk,3),iz(nk,3)
C
      WRITE(*,*) 'Enter n,m for C_h.   e.g., 10 5'
      READ(*,*) n,m
      WRITE(*,*) 'Enter the number of the unit cell N_u = '
      READ(*,*) nu
      CALL UNITCEL(n,m,np,nq,ndr)
C
C     OPEN(61,file='unitube.d')
C     READ(61,*) nn
C     READ(61,*) t
C     DO 10 i=1,nn
C        READ(61,*) j,x(i),y(i),z(i)
C  10 CONTINUE
C     CLOSE(61,status='keep')
C
      WRITE(*,*) 'Number of total atoms = ',nu*nn
      IF(nu*nn.GT.nk) Stop 'Change parameter nk to a larger value.'
C
      ii=nn
      DO 20 i=1,nn
         DO 21 jj=0,nu-1
            iii=ii+i+jj*ii
            z(iii)=z(i)+dfloat(jj+1)*t
            x(iii)=x(i)
            y(iii)=y(i)
   21    CONTINUE
   20 CONTINUE
C
C     OPEN(60,FILE='cntube.xyz')
C     WRITE(60,*) nu*nn
C     WRITE(60,*) ' '
C     DO 100 i=1,nu*nn
C            WRITE (60,700) x(i),y(i),z(i)
C 100 CONTINUE
  700 FORMAT('C ',3F10.5)
C     CLOSE(60)
C
      STOP
      END

C **********************************************************************

      SUBROUTINE GEN11(n,m,np,nq,ndr)
C
C  Parameters for a (n,m) nanotube
C
C  Made by 10/16/95 R. Saito
C
C  Reference: "Physical Properties of Carbon Nanotubes"
C              Imperial College Press,
C              by R. Saito, G. Dresselhaus and M. S. Dresselhaus
C
C
      dimension nnp(1000), nnq(1000)
C
      itest=1
      itest1=0
C
C      WRITE(*,*) 'Enter n,m for C_h.   e.g., 10 5'
C      READ(*,*) n,m
      nd=igcm(n,m)
      IF(mod((n-m),3*nd).EQ.0) THEN
      ndr=3*nd
      ELSE
      ndr=nd
      ENDIF
      IF(itest.EQ.1) WRITE(*,21) n,m
   21 FORMAT(' C_h = (',I3,',',I3,')','     : CHIRAL VECTOR')
      IF(itest.EQ.1) WRITE(*,22) nd
   22 FORMAT('   d =  ',I3)
      IF(itest.EQ.1) WRITE(*,23) ndr
   23 FORMAT(' d_R =  ',I3)
C
      a=sqrt(3.0)*1.42
      eps=1.0E-5
C  L
      l2=n*n+m*m+n*m
      If(l2.LE.0) stop 'l2.LE.0'
      l=int(sqrt(dfloat(l2))+eps)
      IF((l2-l**2).EQ.0) THEN
      IF(itest.EQ.1) WRITE(*,31) l
   31 FORMAT(' L/a = ',I4)
      ELSE
      IF(itest.EQ.1) WRITE(*,32) l2
   32 FORMAT(' L/a = square root of',I6)
      ENDIF
      dt=a*sqrt(dfloat(l2))/3.1415926525
      rt= dt*0.5
      IF(itest.EQ.1) WRITE(*,33) dt,rt
   33 FORMAT(' d_t = ',F10.5,' A    r_t = ',F10.5,' A   : DIAMETER ', 
     & '& RADIUS')
C  T
      nr=(2*m+n)/ndr
      ns=-(2*n+m)/ndr
      IF(itest.EQ.1) WRITE(*,41) nr,ns
   41 FORMAT('   T = (',I3,',',I3,')')
      nt2=3*l2/ndr/ndr
      nt=int(sqrt(dfloat(nt2))+eps)
      IF((nt2-nt**2).EQ.0) THEN
      IF(itest.EQ.1) WRITE(*,42) nt
   42 FORMAT(' T/a = ',I3)
      ELSE
      IF(itest.EQ.1) WRITE(*,43) nt2
   43 FORMAT(' T/a = square root of',I6)
      ENDIF
C  N
      nn=2*l2/ndr
      IF(itest.EQ.1) WRITE(*,51) nn
   51 FORMAT('   N = ',I4)
C  R
      ichk=0
      IF(nr.EQ.0) THEN
         n60=1
         ELSE
         n60=nr
      ENDIF
C
C      itest2=0
      itest2=1
C
      DO 61 np=-abs(n60),abs(n60)
         DO 62 nq=-abs(ns),abs(ns)
               j2=nr*nq-ns*np
               IF(j2.EQ.1) THEN
               j1=m*np-n*nq
               IF(itest2.EQ.1) WRITE(*,68) n,m,nr,ns,np,nq,j1,j2
   68          FORMAT('n, m, nr, ns, np, nq, j1, j2 = ',6I4,2I6)
               IF((j1.GT.0).AND.(j1.LT.nn)) THEN
               ichk=ichk+1
               nnp(ichk)=np
               nnq(ichk)=nq
               ENDIF
               ENDIF
   62    CONTINUE
   61 CONTINUE
C
      IF(ichk.EQ.0) THEN
      WRITE(*,*) 'n, m, nr, ns = ', n,m,nr,ns
      STOP 'not found p, q strange!!'
      ENDIF
C
      itest3=1
C
      IF(ichk.GE.2) THEN
      IF(itest3.EQ.1) THEN
      WRITE(*,*) 'n, m, nr, ns = ',n,m,nr,ns
      WRITE(*,*) 'ichk = ',ichk,' ndr = ',ndr
      WRITE(*,67) (nnp(i),nnq(i),(m*nnp(i)-n*nnq(i)),nn,i=1,ichk)
      STOP 'more than 1 pair of p, q strange!!'
      ENDIF
C
      IF((nr.NE.0).AND.(ns.NE.0)) THEN
C
      IF(itest1.EQ.1) THEN
      DO 77 i=1,ichk
      IF((m*nnp(i)-n*nnq(i)).LT.nn) GOTO 777
   77 CONTINUE
      ENDIF
C
      WRITE(*,*) 'n, m, nr, ns = ',n,m,nr,ns
      WRITE(*,*) 'ichk = ',ichk,' ndr = ',ndr
      WRITE(*,67) (nnp(i),nnq(i),(m*nnp(i)-n*nnq(i)),nn,i=1,ichk)
   67 FORMAT(' (',I2,',',I2,')  mp-nq = ',I3,'  N = ',I3)
      STOP ' more than 1 pair of p, q strange!!'
      ENDIF
      ENDIF
 777  CONTINUE
C
      IF(itest.EQ.1) WRITE(*,63) nnp(1),nnq(1)
   63 FORMAT('  R = (',I3,',',I3,')')
C  M = mp - nq
      mmm=m*nnp(1)-n*nnq(1)
      WRITE(*,69) mmm
   69 FORMAT('  M = ',I3)
C 
      IF(itest.EQ.1) THEN
      np=nnp(1)
      nq=nnq(1)
      ENDIF
    2 CONTINUE
      RETURN
    1 CONTINUE 
      STOP
      END

C
C  calculate the highest common divisor
C
      integer function igcm(ii,jj)
      i=abs(ii)
      j=abs(jj)
      IF(j.GT.1) THEN
      iw=j
      j=i
      i=iw
      ENDIF
      IF(j.EQ.0) THEN
      igcm=i
      RETURN
      ENDIF
      
   10 CONTINUE
      ir=mod(i,j)
      IF(ir.EQ.0) THEN
      igcm=j
      RETURN
      ELSE
      i=j
      j=ir
      ENDIF
      GOTO 10
      END

C **********************************************************************
      SUBROUTINE UNITCEL(n,m,np,nq,ndr)
C
C  Atomic coordiates in the unit cell of a (n,m) nanotube
C
C  Input: n, m
C  Output: unitube.xyz (Coordinates of unitcell)
C          unitube.d (a work file for program mk_cnt)
      IMPLICIT REAL*8(a-h,o-z)
C
C   acc = C-C bond distance
C
      PARAMETER(acc=1.42d0)
      PARAMETER(nk=20000)
      DIMENSION x(nk),y(nk),z(nk)
C
C      WRITE(*,*) 'Enter n, m for C_h=(n,m)'
C      READ(*,*) n,m
C

       CALL gen11(n,m,np,nq,ndr)
C       WRITE(*,*) ' R = (',np,',',nq,')'
C  a = length of unit vector
C  pi = 3.141592    sq3=1.732
      sq3=sqrt(3.0d0)
      pi=4.0d0*atan(1.0d0)
      a=sqrt(3.0d0)*acc
C
C  r=|R|, c=|C_h|, t=|T|
C
      r=a*sqrt(dfloat(np*np+nq*nq+np*nq))
      c=a*sqrt(dfloat(n*n+m*m+n*m))
      t=sqrt(3.0d0)*c/dfloat(ndr)
      WRITE(*,*) '  t = ',t, dfloat(n*n+m*m+n*m)
C      WRITE(*,190) c/2.0d0/pi
C 190 FORMAT(' nradius = ',F10.4)
C
C  nn: the number of hexagon in the unit cell N
      nn = 2*(n**2+m**2+n*m)/ndr
C
C     WRITE(*,*) '   N = ',nn*2
C
C  rs: radius of tube
C
      IF(2*nn.GT.nk) stop 'parameter nk is too small.'
C
      rs=c/(2.0d0*pi)
C
C  q1: the chiral angle for C_h
C  q2: the 'chiral' angle for R
C  q3: the angle between C_h and R
C
      q1=atan((sq3*dfloat(m))/dfloat(2*n+m))
      q2=atan((sq3*dfloat(nq))/dfloat(2*np+nq))
      q3=q1-q2
C
C  q4: a period of an angle for the A atom
C  q5: the difference of the angle between the A and B atoms
C
      q4=2.0d0*pi/dfloat(nn)
      q5=acc*cos((pi/6.0d0)-q1)/c*2.0d0*pi
C
C  h1:
C  h2: Delta z between the A and B atoms
C
      h1=abs(t)/abs(sin(q3))
      h2=acc*sin((pi/6.0d0)-q1)
C
      WRITE(*,*) 'q1: =',q1*180.0d0/pi,'  : CHIRAL ANGLE'
C      WRITE(*,*) 'q2: =',q2*180.0d0/pi
C      WRITE(*,*) 'q4: =',q4*180.0d0/pi
C      WRITE(*,*) 'q5: =',q5*180.0d0/pi
C
C   calculate 2*nn atoms in the unit cell.
C
C The A atom
C
      ii=0
      DO 100 i=0,nn-1
      x1=0
      y1=0
      z1=0
      k=int(dfloat(i)*abs(r)/h1)
      x1=rs*cos(dfloat(i)*q4)
      y1=rs*sin(dfloat(i)*q4)
      z1=(dfloat(i)*abs(r)-dfloat(k)*h1)*sin(q3)
C      z1=dfloat((dfloat(i)*abs(r)-dfloat(k)*h1))*sin(q3)
      kk2=abs(int(z1/t))+1
C
C  Check the A atom is in the unit cell 0 < z1 < t.
C
      IF(z1.GT.t-0.02) THEN
         z1=z1-t*dfloat(kk2)
      ENDIF
      IF(z1.LT.-0.02) THEN
         z1=z1+t*dfloat(kk2)
      ENDIF
      ii=ii+1
      x(ii)=x1
      y(ii)=y1
      z(ii)=z1
C      WRITE(*,*) "z(ii) = ",z(ii)
C
C  The B atoms
C
      z3=(dfloat(i)*abs(r)-dfloat(k)*h1)*sin(q3)-h2
      ii=ii+1
C
C  Check the B atom is in the unit cell 0 < z3 < t.
C
      IF((z3.GE.-0.02).AND.(z3.LE.t-0.02)) THEN
C yes
         x2=rs*cos(dfloat(i)*q4+q5)
         y2=rs*sin(dfloat(i)*q4+q5)
C         z2=dfloat(dfloat(i)*abs(r)-dfloat(k)*h1)*sin(q3)-h2
         z2=(dfloat(i)*abs(r)-dfloat(k)*h1)*sin(q3)-h2
         x(ii)=x2
         y(ii)=y2
         z(ii)=z2
      ELSE
C no
         x2=rs*cos(dfloat(i)*q4+q5)
         y2=rs*sin(dfloat(i)*q4+q5)
C         z2=dfloat(dfloat(i)*abs(r)-dfloat(k+1)*h1)*sin(q3)-h2
         z2=(dfloat(i)*abs(r)-dfloat(k+1)*h1)*sin(q3)-h2
         kk=abs(int(z2/t))+1
      IF(z2.GT.t-0.02) THEN
         z2=z2-t*dfloat(kk)
      ENDIF
      IF(z2.LT.-0.02) THEN
         z2=z2+t*dfloat(kk)
      ENDIF
         x(ii)=x2
         y(ii)=y2
         z(ii)=z2
      ENDIF
  100 CONTINUE
C
C  Output to the file tube.xyz
C
C     OPEN(60,file='unitube.xyz')
C     WRITE(60,*) 2*nn
C     WRITE(60,*) ' '
C     DO i=1,nn*2
C        WRITE(60,117) x(i),y(i),z(i)
C     ENDDO
C 117 FORMAT('C',3F10.5)
C     CLOSE(60)
C
C Output to the file tube.xyz
C
C     OPEN(60,file='unitube.d')
C     WRITE(60,*) 2*nn
C     WRITE(60,118) t,acc
C 118 FORMAT(2F25.5)
C     DO i=1,nn*2
C        WRITE(60,116) i,x(i),y(i),z(i)
C     ENDDO
C 116 FORMAT(I5,3F25.20)
C     CLOSE(60)

C
C Output to the fdf file 
C
C     OPEN(60,file='nanotube.fdf')
C     WRITE(60,*)""
C
C     write(60,*) 'SystemLabel   Nanotube'
C     write(60,*) 'NumberOfSpecies    1'
C     write(60,'(A14,I4)') 'NumberOfAtoms', 2*nn
C     write(60,*) 'SpeciesAndOrbitals'
C     write(60,*) 1, 1
C     write(60,*) 'FileParameters     c-c-pz'
C     write(60,*) '%block ChemicalSpeciesLabel'
C     write(60,*) '  1  6  C'
C     write(60,*) '%endblock ChemicalSpeciesLabel'
C     write(60,*) 'NumberOfLines   2'
C     write(60,'(A14,I5)') 'NumberOfBands', 2*nn
C     write(60,'(A21,I5)') 'NumberOfOccupiedBands', nn
C     write(60,*) 'NumberOfDivisions 50'
C     write(60,*) 'KPointsScale ReciprocalLatticeVectors  '
C     write(60,*) 'KPointsAndLabels  '                                
C     write(60,*) ' 0.0000000       0.0000000       0.5000000 X '
C     write(60,*) ' 0.0000000       0.0000000       0.0000000 G '    
C     write(60,*) ' 0.0000000       0.0000000       0.5000000 X '    
C     write(60,*) 'SlaterKoster'                                      
C     write(60,*) 'MaximumDistance 1.5'                               
C     write(60,*) 'InputEnergiesInEV  '                               
C     write(60,*) 'LatticeConstant       1.0 Ang'
C     write(60,*) '%block LatticeVectors'
C     write(60,'(3F10.5)')   rs+20.0, 0.d0 ,0.d0
C     write(60,'(3F10.5)') 0.d0, rs+20.0, 0.d0
C     write(60,'(3F10.5)') 0.d0, 0.d0, t
C     write(60,*) '%endblock LatticeVectors'
C     write(60,*) 'AtomicCoordinatesFormat   ScaledCartesianAng'
C     write(60,*) '%block AtomicCoordinatesAndAtomicSpecies' 
C     DO i=1,nn*2
C        WRITE(60,120) x(i),y(i),z(i), 1
C     ENDDO
C     write(60,*) '%endblock AtomicCoordinatesAndAtomicSpecies' 
C120  FORMAT(3F17.9,I5)
C     CLOSE(60)

      OPEN(60,file='nanotube.cif')
      WRITE(60,*) 'data_block_1'
      WRITE(60,121) '_cell_length_a', t 
      WRITE(60,121) '_cell_length_b', rs+20.0
      WRITE(60,121) '_cell_length_c', rs+20.0
      WRITE(60,121) '_cell_angle_alpha', 90.0
      WRITE(60,121) '_cell_angle_beta', 90.0
      WRITE(60,121) '_cell_angle_gamma', 90.0
      WRITE(60,*) 'loop_'
      WRITE(60,*) '_atom_site_type_symbol'
      WRITE(60,*) '_atom_site_fract_x'
      WRITE(60,*) '_atom_site_fract_y'
      WRITE(60,*) '_atom_site_fract_z'
  121 FORMAT(A17, F25.5)
      DO i=1,nn*2
         WRITE(60,122) 'C',z(i)/t,(x(i)+10+rs/2.0)/(rs+20.0), 
     .        (y(i)+10+rs/2)/(rs+20)
      ENDDO
  122 FORMAT(A2, 3F17.10)
      CLOSE(60)

      RETURN
      STOP
      END


