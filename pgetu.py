import numpy as np
import scipy;from scipy import ndimage
import matplotlib;from matplotlib import pyplot as plt

GX,GY,GB,GE,GPA,GSM,GSA = range(7)
SX,SY,SF,SW,SE,SPA = range(6)

def shift_func(coords,xoff,yoff):
    return (yoff[coords[0],coords[1]],xoff[coords[0],coords[1]])

def sie (x,y,lb,le,lpad,smag,sangd):
    lpa,sang = np.deg2rad(lpad),np.deg2rad(sangd)
    cpa,spa = np.cos(lpa),np.sin(lpa)
    xp,yp =  x*cpa + y*spa, -x*spa + y*cpa
    phis = np.arctan2 (yp,xp)
    sphis,cphis = np.sin(phis),np.cos(phis)
    axrat = 1.0-le
    sqaxrat = np.sqrt(axrat)
    axprime = np.sqrt(1.0-axrat*axrat)
    kappa = 0.5*sqaxrat/np.hypot(xp/lb,axrat*yp/lb)
    gamma1 = kappa*(2.*sphis*sphis-1.0)
    gamma2 = -kappa*2.*sphis*cphis
    fasin = np.arcsin(axprime*sphis)
    fasinh = np.arcsinh (axprime*cphis/axrat)
    fratio = sqaxrat/axprime
    pot = np.hypot(x,y)*lb*fratio*(sphis*fasin+cphis*fasinh)
    xs,ys = xp - lb*fratio*fasinh, yp - lb*fratio*fasin
    xp,yp = xs*cpa - ys*spa, xs*spa + ys*cpa
    csang,ssang = np.cos(2.*sang),np.sin(2.*sang)
    csath,ssath = np.cos(2.*(sang-lpa)),np.sin(2.*(sang-lpa))
    xp += -smag * (x*csang + y*ssang)
    yp += -smag * (x*ssang - y*csang)
    gamma1 += smag*csath
    gamma2 -= smag*ssath
    return xp,yp,kappa,gamma1,gamma2,pot

def lens (g,a=[],gpix=128):
    try:
        x,y = a[1]-g[GX],a[0]-g[GY]
    except:
        x,y=np.meshgrid(np.arange(gpix)-g[GX],np.arange(gpix)-g[GY]) 
    xp,yp,kappa,gamma1,gamma2,pot = sie (x,y,g[GB],\
                              max(g[GE],1.E-5),g[GPA],g[GSM],g[GSA])
    u = np.array(([kappa+gamma1, gamma2, gamma2, kappa-gamma1]))
    return u,np.array([x-xp,y-yp]),pot

def mxrot (a, angle):
    t = np.zeros_like(a)
    c,s = np.cos(np.deg2rad(angle)),np.sin(np.deg2rad(angle))
    t[0],t[1],t[2],t[3] = c*a[0]+s*a[2], c*a[1]+s*a[3],\
                         -s*a[0]+c*a[2],-s*a[1]+c*a[3]
    return t

def scene (ag, s, gpix):
    ag = [ag] if ag.ndim==1 else ag
    aoff = np.zeros((2,s.shape[0],s.shape[1]))
    for g in ag:
        u,alpha,pot = lens(g,[],gpix)
        aoff[0] += alpha[0]
        aoff[1] += alpha[1]
    xoff,yoff = np.meshgrid (np.arange(s.shape[0],dtype='float'),\
                             np.arange(s.shape[1],dtype='float'))
    xoff -= aoff[0]
    yoff -= aoff[1]
    return ndimage.geometric_transform(s,shift_func,extra_arguments=(xoff,yoff))

def gcsm (ag,src,gpix):
    if len(src):
        allpot,a,aoff = 0.0, np.array([1.,0.,0.,1.]), np.array([0.0,0.0])
    else:
        allpot = np.zeros((gpix,gpix))
        a = np.zeros((4,gpix,gpix)); a[0]=1.; a[3]=1.
        aoff = np.zeros((2,gpix,gpix))

    for g in ag:
        u,alpha,pot = lens(g,[src[0],src[1]],gpix) if len(src) else lens(g,[],gpix)
        allpot += pot
        kappa = 0.5*(u[0]+u[3])
        u[0]-=kappa
        u[3]-=kappa
        u = mxrot(u,-2*g[GPA])
        u[0]+=kappa
        u[3]+=kappa
        a -= u
        aoff[0] += alpha[0]
        aoff[1] += alpha[1]
#    tdel -= allpot           # !!! tdel comes out wrong at the moment
    tdel = 0.0
    mag = 1.0/(a[0]*a[3]-a[1]*a[2])
    return a,aoff,allpot,tdel

def draw_caus (s,ag,filt=2.0,blc=[0,0],trc=[0,0],gpix=128):
    
    trc =  [gpix,gpix]#trc if trc.sum() else
    a,aoff,allpot,tdel = gcsm (ag,[],gpix)
    mag = (a[0]*a[3]-a[1]*a[2])
  #  plt.clf()
  #  plt.imshow(mag)
  #  plt.show()
 #   plt.contour(scipy.ndimage.gaussian_filter(mag[blc[1]:trc[1],blc[0]:trc[0]],\
     #             filt),levels=[0.0],colors='black')
    print mag.shape
    copymag = np.copy(mag)
    zoomfac = 3
  #  copymag = scipy.ndimage.zoom(copymag, zoomfac) 
    print mag.shape
    caus = 0.0*mag
    print caus.shape
    
    for iy in range(gpix):
        for ix in range(gpix):
            if mag[iy,ix]**2. < 0.0001:
                jy = iy-aoff[1,iy,ix]
                jx = ix-aoff[0,iy,ix]
                try:
                    caus[np.int(jy),np.int(jx)]=1.0
                except:
                    pass
  #  caus = scipy.ndimage.zoom(caus, zoomfac)
#    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    trc[0]*=zoomfac
    trc[1]*=zoomfac
    print trc, blc
    magmag = 1/np.abs(mag)
    np.putmask(magmag, magmag>10000, 10000)
  #  magmag = np.log10(magmag)
 #   ax.imshow(magmag)
   # ax.imshow(caus[blc[1]:trc[1],blc[0]:trc[0]], cmap = 'gray_r')
 #   ax.contour(scipy.ndimage.gaussian_filter(copymag[blc[1]:trc[1],blc[0]:trc[0]],\
 #                 filt),levels=[0.0],colors='black')
 #   from matplotlib.patches import Ellipse              
 #   ax.add_patch(Ellipse((s[0],s[1]), width=1, height = 1, fill = 0, color = 'red'))
 #   for j in range(1):
 #       for i in xrange(comps[j].shape[0]):
 #           xpos = comps[j][i,0]
 #           ypos = comps[j][i,1]
 #           ax.add_patch(Ellipse((xpos,ypos), width=5, height = 5, fill = 0, color = 'blue'))
 #   plt.show()         
    return caus, magmag, copymag  , blc, trc 
    

def groot(ag,s):
    return s
