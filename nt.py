import numpy as np
from math import gcd
import supertb as tb


def gen11(n,m,cc=1.42):
    
    nd = int(gcd(n,m))
    if int(n-m)%int(3*nd) == 0:
        ndr=3*nd
    else:
        ndr=nd
        
    #print ('nd ', nd, ndr)
        
    a = np.sqrt(3.)*cc
    eps = 1.0e-5
    l2=n*n+m*m+n*m
    l=int(np.sqrt(float(l2))+eps)
    #if l2-l**2==0:
    #    print ('L/a = ', l)
    #else:
    #    print ('L/a = square root of', l2)
    dt=a*np.sqrt(float(l2))/3.1415926525
    rt= dt*0.5
    #print ('dt = ', dt, ' rt = ', rt)
    
    nr = (2*m+n)/ndr
    ns = -(2*n+m)/ndr
    nn=2*l2/ndr
    #print ('T = ', nr, ns)
    #print ('N =', nn)
    
    ichk = 0
    if nr==0:
        n60=1
    else:
        n60=nr
    
    nnp = []
    nnq = []
    for npi in range(-int(np.abs(n60)),int(np.abs(n60))+1):
        for nqi in range(-int(np.abs(ns)),int(np.abs(ns))+1):
            j2 = nr*nqi-ns*npi
            if (j2==1):
                j1=m*npi-n*nqi
                #print (n,m,nr,ns,npi,nqi,j1,j2)
                if (j1>0) and (j1<nn):
                    nnp.append(npi)
                    nnq.append(nqi)
                    ichk+=1
    #print (ichk)         
    npi = nnp[0]
    nqi = nnq[0]
    
    #print ('R =', npi, nqi)
    return npi, nqi, ndr
    
def create_nt(n,m,acc=1.42,nk=1000):
    
    npi, nqi, ndr = gen11(n,m,cc=acc)
    #print (npi,nqi,ndr)
    sq3 = np.sqrt(3.)
    a = sq3*acc
    
    r = a*np.sqrt(float(npi*npi+nqi*nqi+npi*nqi))
    c = a*np.sqrt(float(n*n+m*m+n*m))
    t = sq3*c/float(ndr)
    #print ('  t = ',t, float(n*n+m*m+n*m))
    
    nn = 2*(n**2+m**2+n*m)/ndr
    
    if (2*nn > nk):
        raise IOError("parameter nk is too small")
        
    rs = c/(2*np.pi)
    q1 = np.arctan(sq3*m/(2*n+m))
    q2 = np.arctan(sq3*nqi/(2*npi+nqi))
    q3 = q1 - q2
    q4 = 2*np.pi/nn
    q5 = acc*np.cos((np.pi/6)-q1)/c*2*np.pi
    
    h1 = np.abs(t)/np.abs(np.sin(q3))
    h2 = acc*np.sin((np.pi/6)-q1)
    
    #print ('q1: =',q1*180.0/np.pi,'  : CHIRAL ANGLE')
    #print ('q2: =',q2*180.0/np.pi,'  : CHIRAL ANGLE')
    #print ('q4: =',q4*180.0/np.pi,'  : CHIRAL ANGLE')
    #print ('q5: =',q5*180.0/np.pi,'  : CHIRAL ANGLE')
    
    ii = 0
    x = []
    y = []
    z = []
    for i in range(int(nn)):
        k = int(i*np.abs(r)/h1)
        x1 = rs*np.cos(i*q4)
        y1 = rs*np.sin(i*q4)
        z1 = (i*np.abs(r)-k*h1)*np.sin(q3)
        kk2 = np.abs(int(z1/t))+1
        
        if (z1 > t-0.02):
            z1 = z1-t*float(kk2)
        if (z1 < -0.02):
            z1 = z1+t*float(kk2)
        x.append(x1)
        y.append(y1)
        z.append(z1)
        
        z3 = (i*np.abs(r)-k*h1)*np.sin(q3)-h2
        if (z3>=-0.02) and (z3<=t-0.02):
            x2 = rs*np.cos(i*q4+q5)
            y2 = rs*np.sin(i*q4+q5)
            z2 = (i*np.abs(r)-k*h1)*np.sin(q3)-h2
        else:
            x2 = rs*np.cos(i*q4+q5)
            y2 = rs*np.sin(i*q4+q5)
            z2 = (i*np.abs(r)-(k+1)*h1)*np.sin(q3)-h2
            kk = np.abs(int(z2/t))+1
            if (z2>t-0.02):
                z2 = z2-t*float(kk)
            if (z2<-0.02):
                z2 = z2+t*float(kk)
        x.append(x2)
        y.append(y2)
        z.append(z2)
        
    #print (len(x))
    #print (nn)
        
    lattice = tb.Lattice.from_parameters(t,rs+20.,rs+20.,90.,90.,90.)
    fracs = np.zeros((2*int(nn),3))
    #carts = np.zeros((2*int(nn),3))
    for i in range(2*int(nn)):
        fracs[i] = [z[i]/t,(x[i]+10+rs/2.)/(rs+20.),(y[i]+10+rs/2)/(rs+20)]
        #carts[i] = [z[i],x[i],y[i]]
    return tb.Structure(lattice,['C']*(2*int(nn)),fracs,coords_are_cartesian=False)
    #return tb.Structure(lattice,['C']*(2*int(nn)),carts,coords_are_cartesian=True)
    
