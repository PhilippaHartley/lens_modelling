#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 11:30:13 2018

@author: phartley
"""
import numpy as np, scipy as sp,os,sys,matplotlib,pgetu,warnings

from scipy import ndimage,optimize; from scipy.optimize import fmin, basinhopping
from scipy.optimize import fmin_l_bfgs_b as lbfgs
from scipy.optimize import fmin_bfgs as bfgs
from scipy.optimize import fmin_powell as powl
import time, pickle, njj, matplotlib.pyplot as plt
from astropy.io.fits import getdata,getheader
from matplotlib import gridspec
from numpy import random
from numpy import log, exp, pi
import scipy.stats, scipy
import pymultinest
import json
#import sep
import time
import sep

GX,GY,GB,GE,GPA,GSM,GSA = range(7); SX,SY,SF,SW,SE,SPA,SSER = range(7)
warnings.simplefilter("ignore", RuntimeWarning)
plt.rcParams['figure.subplot.wspace']=plt.rcParams['figure.subplot.hspace']=0.01
plt.rcParams['image.origin']='lower'
plt.rcParams['image.interpolation']='nearest'
figs='.png'




def extract (a_stack, comps,x0, lin,ospar,ogpar):
    timer = time.time()
    npix = 1200/2
    from matplotlib.patches import Ellipse
    global new_plotname
    global last_retval
    fig = plt.figure(1, figsize = (12, 7))
    ax = plt.subplot(121)
    both = np.zeros((npix,npix))
    if a_stack.ndim == 3:
        dim = 2
        for j in range(dim):
            both += a_stack[:,:,j]
    else:
        both  = a_stack
        dim = 1
    plt.imshow(both, cmap = 'gray_r')
    retval = 0
    ax = plt.subplot(122)
    imsize = 1200/2
    im = np.zeros((imsize,imsize))
    implot = plt.imshow(im, cmap = 'gray_r')
    for j in range(dim):
        if j ==0:
            col = '#c63971'
        else: 
            col = '#5c7ad6'
        if dim == 2:
            a = a_stack[:,:,j]
        else:
            a = a_stack  
        noise = random.normal(0.,0.0000000001,(npix,npix)) # needed to make sextractor behave
        a+=noise
        newtime = time.time()
        timer = np.copy(newtime)
        a = a.copy(order = 'C')
        thresh = 1.5 * 0.0000000001
        objects = sep.extract(a, thresh)
        len(objects)  # number of objects
        cat = np.array([objects['npix'], objects['flux'],objects['thresh'], objects['x'],objects['y']]).T
        if len(cat)==0:
            plt.clf()
            return 2e6, last_retval
        if cat.ndim==1:
            cat = np.array([cat])
            cat = cat[cat[:,1]>1e-6]
        if cat.shape[0] <comps[j].shape[0]:
            plt.clf()
            return 1e6, last_retval
        sortcat = cat[np.argsort(cat[:,1]),:][::-1]
        plt.yticks([])
        elsi = 10
        for i in xrange(comps[j].shape[0]):
            xpos = comps[j][i,0]
            ypos = comps[j][i,1]
            ax.add_patch(Ellipse((xpos,ypos), width=elsi, height = elsi, fill = 0, color = col))
        for i in xrange(comps[j].shape[0]):
            xpos = sortcat[i,3]
            ypos = sortcat[i,4]
            flux = sortcat[i,1]*4000
            ax.add_patch(Ellipse((xpos,ypos), width=elsi, height = elsi, fill = 1, color = col))
            model_comp = np.array([xpos,ypos])
            try:
                model_comps = np.vstack((model_comps, model_comp))
            except:
                model_comps = model_comp
        for i in comps[j]:
            i_tiled =  np.tile(i, comps[j].shape[0]).reshape(comps[j].shape[0],2)
            legs = np.vstack((i_tiled[:,0]-model_comps[:,0],i_tiled[:,1]-model_comps[:,1])).T
            distances = np.hypot(legs[:,0],legs[:,1])
            try:
                all_dist = np.vstack((all_dist,distances))
            except:
                all_dist = distances
        for k in range(len(all_dist)):
            if k == 3:            
                err = 2.
            else: 
                err = 4.    
            retval += np.min(all_dist[k])**2/err
            delarg = np.argwhere(all_dist[k] == np.min(all_dist[k]))
            all_dist = np.delete(all_dist,delarg,1)
         #   print 'retval: ', retval             
        del all_dist
        del model_comps      
    ax.set_title(retval)
    fig.text(0.,0.01,x0)
    fig.text(0.,0.85,lin)
    new_dirname = 'sextractor_plots'
    pltname='/best_fit'
    #plt.show()
    if not last_retval == None:
        if retval < last_retval:
            last_retval = np.copy(retval)    
            print ('retval', retval)
            print( 'plotting')
            plt.yticks([])
            plt.savefig(new_plotname)
            plt.clf()
            os.system('rm par.log')
            f=open('par.log','a')
            count = 0
            np.save('ospar', ospar)
            np.save('ogpar', ogpar)
            for i in ospar:
                if len(ospar)==1:
                    f.write( '    spar1 = np.array([[%f, %f, %f, %f, %f, %f]])' %(i[0],i[1],i[2],i[3],i[4],i[5])+'\n')
                elif len(ospar)==2:
                    if count==0:
                        f.write( '    spar1 = np.array([[%f, %f, %f, %f, %f, %f],' %(i[0],i[1],i[2],i[3],i[4],i[5]))
                    if count==1:
                        f.write( '[%f, %f, %f, %f, %f, %f]])' %(i[0],i[1],i[2],i[3],i[4],i[5])+'\n')
                count+=1            
            for i in ogpar:
                f.write( '    gpar1 = np.array([[%f , %f , %f , %f , %f , %f, %f]])' %(i[0],i[1],i[2],i[3],i[4],i[5],i[6]))
            f.close()
        else:
            plt.clf()
    else:
        n = 0
        while os.path.exists('{}{:d}.png'.format(new_dirname+pltname, n)):
            n += 1
        new_plotname = '{}{:d}.png'.format(new_dirname+pltname, n)
        last_retval = 2.0e6
    return    retval, last_retval
        
def mksersic (dist, flux, re, axrat=1.0, angle=0.0, sser=1.0):
    sinth = np.sin(angle*np.pi/180.0)
    costh = np.cos(angle*np.pi/180.0)
    r = np.array([-sinth,costh,-costh,-sinth])
    rt = np.array([-sinth,-costh,costh,-sinth])
    sig = np.array([re,0.0,0.0,re*axrat])
    scr1 = mxmul (sig,r)
    scr2 = mxmul (rt, scr1)
    scr1 = mxinv (scr2)
    ex = scr1[0]*dist[0]+scr1[1]*dist[1]
    ey = scr1[2]*dist[0]+scr1[3]*dist[1]
    return (flux/axrat)*np.exp(-(pow((np.hypot(ex,ey)/re),1./sser)))

def get_splane (s,cen,ssiz):
    x1,x2,y1,y2 = cen[0]-ssiz/2.,cen[0]+ssiz/2.,cen[1]-ssiz/2.,cen[1]+ssiz/2.
    a=np.zeros((spsiz,spsiz))
    for i in s:
        pos = [spsiz*(i[SX]-x1)/(x2-x1),spsiz*(i[SY]-y1)/(y2-y1)]
        a += njj.mkgauss ([spsiz,spsiz],pos,i[SF],i[SW]*float(spsiz)/float(ssiz),i[SE],i[SPA],10.0)
    return a

def replace_x0 (spar,gpar,sopt,gopt,x0):
    lin = np.append(np.ravel(spar),np.ravel(gpar))
    olin = np.append(np.ravel(sopt),np.ravel(gopt))
    lin[olin] = x0
    osrc = lin[:len(np.ravel(spar))].reshape(spar.shape)
    ogpar = lin[len(np.ravel(spar)):].reshape(gpar.shape)
    return osrc,ogpar

def lensim (g,s,psf, gpix):
    u,alpha,pot,tdel = pgetu.gcsm(g,[],gpix)
   # pgetu.draw_caus (ag,filt=2.0,blc=[0,0],trc=[0,0],gpix=128)
    x,y=np.meshgrid(np.arange(gpix),np.arange(gpix))
    dist=np.array([x-alpha[0]-s[SX],y-alpha[1]-s[SY]])
    if len(s) == 6:
        a=njj.mkpgauss(dist,s[SF],s[SW],s[SE],s[SPA])
    elif len(s) == 7:
        a=mksersic(dist,s[SF],s[SW],s[SE],s[SPA],s[SSER])
    if psf['fwhm'] > 0.0:
        #try a smaller kernel (not size of input image) to avoid oom
        kpix = 50
        psf=njj.mkgauss([kpix,kpix],[kpix/2.,kpix/2.],1.0,psf['fwhm'],psf['axrat'],psf['pa'])
        a=sp.ndimage.convolve(a,psf)
    return a


        
        
        
        
def plotall (spar, a, b, c,x0, retval):


        
        gs = gridspec.GridSpec(1,3,width_ratios=[1.,1.8,1], wspace = 0.1)

        pscale = 1./(3600.*abs(cdelt))
        cen,ss=np.array([spar[0,0],spar[0,1]]),3.0
        imsrc = get_splane(spar,cen,ss)   # ss=image-pix covered by src plane9.02799783e-01,
        fig = plt.figure(1,(15,10))
        cm = 'jet'
        
        gs0 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
        plt.subplot(gs0[0],xticks=[],yticks=[])
        plt.imshow(ndimage.zoom(imsrc*1.0e6,8), cmap = cm)
        plt.colorbar(shrink=0.25)
        plt.rcParams['font.size'] = 12
        sscale = 1.0
       # b = b**2     
        while pscale*sscale/ss > 0.8:
            sscale /= 10.0
       # print 'bbar', bbar , bbar+pscale*sscale*gpix/ss
        plt.plot([bbar,bbar+pscale*sscale*gpix/ss],[bbar,bbar],'w-')
        gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1], hspace=0, wspace = 0,width_ratios=[1,1.25])
        plt.subplot(gs1[0],xticks=[],yticks=[])
        plt.imshow(ndimage.zoom(b*1.0e6,8),cmap = cm,vmin=1.e6*b.min(),vmax=1.e6*b.max())
        plt.plot([bbar,bbar+pscale],[bbar,bbar],'w-')
        plt.subplot(gs1[1],xticks=[],yticks=[])
        plt.imshow(ndimage.zoom(a*1.e6,8),vmin=1.e6*b.min(),vmax=1.e6*b.max(), cmap = cm )
        plt.colorbar( shrink = 0.25)
        plt.plot([bbar,bbar+pscale],[bbar,bbar],'w-')
        gs2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[2], hspace=0, wspace = 0)
        plt.subplot(gs2[0],xticks=[],yticks=[])
        plt.imshow(ndimage.zoom(1.0e6*(b-a),8),vmin=1.e6*b.min(),vmax=1.e5*b.max(), cmap = cm)
        plt.colorbar(shrink=0.25)
        plt.plot([bbar,bbar+pscale],[bbar,bbar],'w-')
        plt.savefig('peluv.pdf',bbox_inches='tight')
        plt.clf()

def optim_func(x0, *x):
    spar,sopt,gpar,gopt,gpix, psfpar,crval,cdelt,crpix,\
                                 b,ismcmc,noise,doplot, comps, lin = x
    psf=dict(); psf['fwhm'],psf['axrat'],psf['pa'] = psfpar
    global sparmin, sparmax, gparmin,gparmin 
    global sparmean, sparstd, gparmean, gparstd
    global last_retval
    ospar,ogpar = replace_x0 (spar,gpar,sopt,gopt,x0)
    for g in ogpar:
        if g[GE]<0.0 or g[GE]>1.0 or g[GSM]<0.0 or g[GSM]>1.0 or g[GB]<0.0:
            return -1.0E6 if ismcmc else 1.0E6
    for i in range(ospar.shape[0]):
        a = lensim (ogpar,ospar[i],psf,gpix)
        try:
            a_stack = np.dstack((a_stack,a))
        except:
            a_stack = a
  #  both= np.zeros((1200,1200))        
  #  for j in range(2):
  #      both += a_stack[:,:,j]
  #  a = both    
        ## now the value to be minimised
   # c = (a-b)**2/noise**2
    usesub = 1
    useextractor =0
    usediff = 0
    if useextractor ==1:
        lin = np.append(np.ravel(spar),np.ravel(gpar))
        retval, last_retval = extract(a_stack, comps,x0, lin,ospar,ogpar)
    if usesub ==1:
        c = (a-b)**2/noise**2
        np.putmask(c,(a<2.*noise)&(b<2.*noise),np.nan)
        pixperbeam = np.pi*psfpar[0]**2*psfpar[1]/(4.*np.log(2.0))
        cpix = len(np.ravel(c)[~np.isnan(np.ravel(c))])
        retval = np.sqrt(np.nansum(c)/cpix) if cpix else 1.0E+6
        if not last_retval == None:
            if retval < last_retval:
                last_retval = retval  
                plotall (ospar, a, b, c,x0, retval)
            else:
                pass
        else:
            last_retval = 1.0e6
            pass        
    ## or use the distance between max pixels in observed and model images
    if usediff ==1: 
        amax = np.argwhere(a==np.max(a))[0]
        bmax = np.argwhere(b==np.max(b))[0]
    #    print( 'locations of max pixels: ',  amax, bmax)
        xdiff = amax[0]-bmax[0]
        ydiff = amax[1]-bmax[1]
        d = np.hypot(xdiff, ydiff)
        retval = d
    penalty = 0.0
    retval = -0.5*(retval+penalty) if ismcmc else retval+penalty    
    f=open('mcmc.log','a')
    for i in x0:
        f.write('%s '%float('%.6g'%i))
    f.write(str(retval)+'\n')
    f.close()
    return retval

def printpar(ospar,ogpar):

    for i in ospar:
        print( '%.4f %.4f %f %.4f %.2f %.2f' % (3600.*cdelt*(i[0]-0.5*gpix),\
              3600.*cdelt*(i[1]-0.5*gpix),\
              i[2],i[3]*3600.*cdelt,i[4],i[5]))
    for i in ogpar:
        print( '%.4f %.4f %f %f %.2f %f %.2f' % (3600.*cdelt*(i[0]-0.5*gpix),\
              3600.*cdelt*(i[1]-0.5*gpix),\
              i[2]*3600.*cdelt,i[3],i[4],i[5],i[6]))
        
        
def model(pos1, width, height1, height2):
    pos2 = pos1 + 0.05
    return  height1 * scipy.stats.norm.pdf(x, pos1, width) + \
        height2 * scipy.stats.norm.pdf(x, pos2, width)

# a more elaborate prior

def prior(cube, ndim, nparams):

    # s params
    cube[0] = 21+cube[0]*20           
    cube[1] = 13+cube[1]*20 # 
    cube[2] = 10**(-8+cube[2]*4)     
    cube[3] = cube[3]*100     
    cube[4] = cube[4]*1    
    cube[5] = cube[5]*180

    #g params
    cube[6] = 23+cube[6]*6      
    cube[7] = 16.+cube[7]*4      
    cube[8] = 5.+cube[8]*20      
    cube[9] =10**(cube[9]*4 - 4)   
    cube[10] = cube[10]*180      
    cube[11] =0+cube[11]*1#10**(cube[11]*4 - 4) 
    cube[12] = cube[12]*180

def loglike(cube, ndim, nparams):
    x0 = []
    for count, i in enumerate(parameters):
        i = cube[count]
        x0 = np.append(x0, i)

    retval = optim_func(x0,spar,sopt,gpar,gopt,gpix, psfpar,crval,cdelt,crpix,\
                      b,False,noise,True, comps, lin)
    loglikelihood = -0.5*(retval)

    return loglikelihood   
    
    
    

def pelim(infile,doopt,spar,sopt,gpar,gopt,gpix, psfpar,\
         crval=[0.,0.],cdelt=[1.,1.],crpix=[0.,0.],noise=1.0, comps=None):
    if len(np.ravel(spar))==12:
        snames = np.array(['SX1','SY1','SF1','SW1','SE1','SPA1']) #,'SX2','SY2','SF2','SW2','SE2','SPA2'
    else: 
        snames = np.array(['SX1','SY1','SF1','SW1','SE1','SPA1'])
    gnames = ['GX','GY','GB','GE','GPA','GSM','GSA']
    global b
    b=getdata(infile)
    
    b=b[0][0] if b.ndim==4 else b
    
    global lin, olin, x0
    
    lin = np.append(np.ravel(spar),np.ravel(gpar))
    olin = np.append(np.ravel(sopt),np.ravel(gopt))
    x0 = lin[olin]
    
    global parameters
    both_names = np.append(np.ravel(snames),np.ravel(gnames))
    parameters = both_names[olin]

    global last_retval
    last_retval = None
    if doopt:
        args=(spar,sopt,gpar,gopt,gpix, psfpar,crval,cdelt,crpix,\
                      b,False,noise,False, comps, lin)
        n_params = len(parameters)
        print ('n_params:', n_params)
        datafile = 'out/1115'
        print (parameters)
        # run MultiNest
        pymultinest.run(loglike, prior, n_params, outputfiles_basename=datafile + '_1_', resume = False, verbose = True)
        json.dump(parameters.tolist(), open(datafile + '_1_params.json', 'w')) # save parameter names
        
        #xopt= lbfgs(optim_func, x0, args=args, approx_grad = 1)
        #xopt= fmin(optim_func, x0, args=args,     maxfun=10000)

    else:   
        optim_func(x0,spar,sopt,gpar,gopt,gpix, psfpar,crval,cdelt,crpix,\
                      b,False,noise,True, comps, lin)
    ospar,ogpar = replace_x0 (spar,gpar,sopt,gopt,xopt if doopt else x0)
#    printpar (ospar,ogpar)

def mcmcinit (x0, nwalkers):
    x0var = np.array([3.,3.,1.e-5,0.02,0.02,1.,1.,1.,0.1,1.e-4,1.0,0.002,1.0])
    p0 = [x0+x0var*(0.5-np.random.rand(len(x0))) for i in xrange(nwalkers)]
    return p0

def mcmc(infile,spar,sopt,gpar,gopt,psfpar,\
         crval=[0.,0.],cdelt=[1.,1.],crpix=[0.,0.],noise=1.0,\
         nwalkers=30, nburn=10, niter=1000):
    b=getdata(infile);  b=b[0][0] if b.ndim==4 else b
    lin = np.append(np.ravel(spar),np.ravel(gpar))
    olin = np.append(np.ravel(sopt),np.ravel(gopt))
    x0 = lin[olin]
    smc = emcee.EnsembleSampler (nwalkers,len(x0),optim_func, args=\
      [spar,sopt,gpar,gopt,psfpar,crval,cdelt,crpix,b,True,noise,False])
    p0 = mcmcinit (x0, nwalkers)
    pos,prob,state = smc.run_mcmc(p0,nburn)
    smc.reset()
    print( '\n\n**** main walking: ****\n\n')
    smc.run_mcmc(pos,niter,rstate0=state)
    return smc
    
if __name__ == "__main__":
    load_params = 0
    factor = 2.
    global last_retval
    
    spar1 = np.array([[ 2.3307e1, 18.250,0.0000892889, 1.307457,7.12575e-2, 1.17286e2]])
    gpar1 = np.array([[  23.1390, 17.50498, 14.8939, 1.257587e-2, 3.85328, 5.0397911e-2, 6.284315e1]])
    
    if load_params == 1:
        sparload = np.load('ospar.npy')
        gparload = np.load('ogpar.npy')
        print ('spar: ', sparload)
        print ('gpar: ', gparload)
        spar1 = sparload
        gpar1 = gparload
    else:
        spar1 = spar1
        gpar1 = gpar1
    
    spar = spar1
    gpar = gpar1        
    
    #   -------- generic parameters for 1115 ------------
    crval,cdelt,crpix=[169.570633  ,  7.766164 ],0.08/3600.,[21,20]
    psfpar=[0.44305/0.08,0.3735/0.44305,-15.45]  # - nb this is in PIXELS: a, axis ratio, angle
   
    y = True
    n = False
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
    gopt,sopt=np.array([y,y,y,y,y,y,y]),\
              np.array([y,y,y,y,y,y])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    
              
    if len(gopt) != 7:
        print ('wrong number of lens booleans!')
        sys.exit()
    if len(sopt)!= 6 and len(sopt) !=12:
        print ('wrong number of source booleans!')
        sys.exit()              
    infile='co_bimage'; gpix=41; doopt=1
    noise=5e-10; spsiz=128; bbar=10
    
    gpar = [gpar] if gpar.ndim==1 else gpar
    spar = [spar] if spar.ndim==1 else spar
    
    if len(np.ravel(sopt)) != len(np.ravel(spar)):
        print ('len(sopt) != len(spar):', len(sopt) ,len(spar) )
        print ('wrong number of source components in par.log!')
        sys.exit()
        
    # This is a list of coordidates of point components
    # It is only needed if the optimisition is being done using the distances between model components and data components (option useextractor in optim_func() ), otherwise pass empty list
    
    #coords:
 #   A1 = [964, 514]
 #   A1a = [939,467]#[944,479]
 #   A2 = [976.1721 ,782.327]
 #   A2a = [1051.,767.]
 #   B1 = [877,353]
 #   B1a=[830.,317.]
 #   B2 = [761,261]
 #   B2a = [725,253]
 #   C1 = [247,214]
 #   C2 = [124,338]
 #   D1 = [379, 1093.1667]
 #   D2 = [323.6432 ,1098.8905]
 #   D2a = [380.6226-56, 1093.1667-27]
 #   D2b = [416,1107]
 #   D2c = [412,1103]
    
    
 #   comps1= np.array([A1a,B2,D2a,C1])/factor#blue
 #   comps2 = np.array([A1,B1,D1,C2])/factor#red
 #   comps = [comps2,comps1]

    comps = []
    
    pelim(infile,doopt,spar,sopt,gpar,gopt,gpix, psfpar,crval,cdelt,crpix,noise,comps)
 
