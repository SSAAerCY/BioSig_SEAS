import numpy as np
import math
import cmath
import scipy.special as special

Pi = math.pi


def mie_abcd(m,x):
    nmax=round(2+x+(4*x**(1./3.)))
    i = 1.0j
    n = np.arange(1,nmax+1,1)
    nu = (n+0.5)
    z = np.multiply(m,x)
    m2 = np.multiply(m,m)
    sqx = np.sqrt(0.5*(math.pi)/x)
    sqz = np.sqrt(0.5*(math.pi)/z)
    bx = (np.multiply((special.jv(nu,x)),sqx))
    bz = np.multiply((special.jv(nu,z)),sqz)
    yx = np.multiply((special.yv(nu,x)),sqx)
    hx = bx+(i*yx)
    sinx = np.array((cmath.sin(x))/x)
    b1x = np.append(sinx,bx[0:int(nmax-1)])
    sinz = np.true_divide(cmath.sin(z),z)
    b1z = np.append(sinz,bz[0:int(nmax-1)])
    cosx = np.array((cmath.cos(x))/x)
    y1x = np.append(-cosx,yx[0:int(nmax-1)])
    h1x = b1x+(i*y1x)
    ax = (np.multiply(x,b1x))-(np.multiply(n,bx))
    az = (np.multiply(z,b1z))-(np.multiply(n,bz))
    ahx = (np.multiply(x,h1x))-(np.multiply(n,hx))
    m2bz = np.multiply(m2,bz)
    antop = (np.multiply(m2bz,ax))-np.multiply(bx,az)
    anbot = (m2*bz*ahx)-(hx*az)
    an = np.true_divide(antop,anbot)
    bn = (bz*ax-bx*az)/(bz*ahx-hx*az)
    cn = (bx*ahx-hx*ax)/(bz*ahx-hx*az)
    dn = m*(bx*ahx-hx*ax)/(m2*bz*ahx-hx*az)
    return np.array([an, bn, cn, dn])

def Mie(m,x):
    if np.any(x)==0: #To avoid a singularity at x=0
        result=0
    else: # This is the normal situation. 
        nmax=round(2+x+(4*x**(1./3.)))
        n1=int(nmax-1);
        n = np.arange(1,nmax+1,1)
        cn=2*n+1
        c1n=np.true_divide((np.multiply(n,(n+2))),(n+1))
        c2n=np.true_divide((np.true_divide(cn,n)),(n+1))
        x2=x*x
        f=mie_abcd(m,x)
        anp=(f[0,:]).real
        anpp=(f[0,:]).imag
        bnp=(f[1,:]).real
        bnpp=(f[1,:]).imag
        g1=np.empty([4,int(nmax)])   
        g1[0,0:int(n1)]=anp[1:int(nmax)] 
        g1[1,0:int(n1)]=anpp[1:int(nmax)]
        g1[2,0:n1]=bnp[1:int(nmax)]
        g1[3,0:n1]=bnpp[1:int(nmax)]
        dn=np.multiply(cn,(anp+bnp))
        q=sum(dn);
        qext=2*q/x2;
        en=np.multiply(cn,(np.multiply(anp,anp)+np.multiply(anpp,anpp)+np.multiply(bnp,bnp)+np.multiply(bnpp,bnpp)))
        q=sum(en);
        qsca=2*q/x2;
        qabs=qext-qsca;
        fn=np.multiply((f[0,:]-f[1,:]),cn)
        gn=(-1)**n;
        f[2,:]=np.multiply(fn,gn)
        q=sum(f[2,:])
        qb=q*(1/q)/x2
        asy1=np.multiply(c1n,(np.multiply(anp,g1[0,:])+np.multiply(anpp,g1[1,:])+np.multiply(bnp,g1[2,:])+np.multiply(bnpp,g1[3,:])))
        asy2=np.multiply(c2n,(np.multiply(anp,bnp)+np.multiply(anpp,bnpp)))
        asy=4/x2*sum(asy1+asy2)/qsca;
        qratio=qb/qsca;
        return qext,qsca,qabs,asy


def q_list(radius,lambdas,index_real,index_imag):
    
    i=1.0j
    reflist = np.array([])
    for real in index_real:
        for j in index_imag:
            realandimag = (real)+(j*i)
            reflist = np.append(reflist,realandimag)
        
    
    xgrid = []
    for lam in lambdas:
        for rad in radius:
            x = np.true_divide((2.0*Pi*rad),float(lam))  
            rl = [rad,lam]
            xgrid.append([x,rad,lam])
            #make a list of xs to deal with later
    
    bigsize = len(xgrid)*len(reflist)
    #keep track of the size needed to calculate everything
    
    result = np.zeros([bigsize,9])
    count = 0
    
    for chunk in xgrid:
        x = chunk[0]
        rad = chunk[1]
        lam = chunk[2]
        for ref in reflist:
            
            qext,qsca,qabs,g = Mie(ref,x)
            sigma_um=(Pi*(rad**2)*(qext))
            #sigma_cm = np.true_divide(sigma_um,10**(8)) 
            result[count,0] = x
            result[count,1] = rad
            result[count,2] = lam
            result[count,3] = ref.real
            result[count,4] = ref.imag
            result[count,5] = qext
            result[count,6] = qsca
            result[count,7] = qabs
            result[count,8] = sigma_um
            count = count+1
    #print 'Size of result = ',len(result)        
    return result[:,8]



