import numpy as np
import scipy.interpolate as interpol
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math as mt
import glob
import sys
import os
import argparse

class fitdata:
    average=np.zeros(1,dtype=np.float64)
    CI=np.zeros(2,dtype=np.float64)
    mode=np.zeros(1,dtype=np.float64)
    imode=np.zeros(1,dtype=np.int)
    sigmamode=np.zeros(2,dtype=np.float64)
    median=np.zeros(1,dtype=np.float64)
    imed=np.zeros(1,dtype=np.int)
    iavg=np.zeros(1,dtype=np.int)
    mean=np.zeros(1,dtype=np.float64)
    imean=np.zeros(1,dtype=np.int)
    sigmamean=np.zeros(2,dtype=np.float64)
    sigmasymmetric=np.zeros(2,dtype=np.float64)
    nmodes=np.zeros(1,dtype=np.int)
    modes=np.zeros(1,dtype=np.float64)
    pdf=np.ndarray([])
    fit=np.ndarray([])
    flags=""

def statsfit(data,datap,cdf,niter,average,smooth=0.0075):
    #create object to contain correct values
    fits=fitdata()
    #first compute cubic spline fit
    splinefit=interpol.splrep(datap,cdf,s=smooth)
    fits.CI=getci(datap,cdf,0.683)
    #now find maximum
    fits.fit=interpol.splev(datap,splinefit)
    firstderivative=interpol.splev(datap,splinefit,der=1)  #stick with old routinefor now, pchip doesn't help,
    secondderivative=interpol.splev(datap,splinefit,der=2) #but at some point need to find a way to constrain the
    mode=datap[np.argmax(firstderivative)]                  #spline such that the first derivative is non-negative
    imode=np.argmax(firstderivative)
    fits.mode=mode
    fits.imode=imode
    fits.pdf=firstderivative
    fits.mean=np.mean(data)
    fits.imean=np.abs(datap-fits.mean).argmin()
    #tbd: find all maxima of curve to flag bi-/multimodal distributions

#    imed=niter/2 #no longer works due to array pruning
    if (niter % 2) == 0: #find median
        imed=np.abs(cdf-0.5).argmin() #find where cdf=0.5 => median
        if cdf[imed] > 0.5: #check that it's not more than half-way through
            imed-=1
        median=(datap[imed]+datap[imed+1])/2.
    else:
        imed=np.abs(cdf-0.5).argmin() #find where cdf=0.5 => median
        if cdf[imed] > 0.5: #check that it's not more than half-way through
            imed-=1
        median=datap[imed]
    fits.median=median
    fits.imed=imed
    if fits.CI[0] == np.min(data):
        fits.flags+='l'
    if fits.CI[1] == np.max(data):
        fits.flags+='u'
    fits.average=fits.mode
    fits.iavg=fits.imode
    if average == 1:
        fits.average=fits.median
        fits.iavg=fits.imed
    if average == 2:
        fits.average=fits.mean
        fits.iavg=fits.imean
    if (average == 0) & (np.isnan(fits.pdf[0])):
        fits.average=fits.median
        fits.iavg=fits.imed
    #return only the important stuff :)
    return fits

def getintervalsym(data,cdf,imed,fraction):
    #for possible future use
    pass

def getinterval(data,cdf,point,fraction):
    #finds the smallest interval of data that contains fraction fraction of the cdf, and point point is inside it.
    interval=np.array([np.min(data),np.max(data)])
    width=interval[1]-interval[0]
    ipoint=np.argmin(np.abs(data-point)) #get location of value that must lie within interval interval
    for i in range(ipoint+1):
        integ=np.float64(0)
        j=np.argmin(np.abs(cdf-cdf[i]-np.float64(fraction)))
        if ((j < len(data)) & (j > ipoint)):
            newwidth=data[j]-data[i]#j-i#
            if newwidth < width:
                width=newwidth
                interval=[data[i],data[j]]

    return interval

def getci(data,cdf,fraction):
    #finds the smallest interval of data that contains fraction fraction of the cdf.
    interval=np.array([np.min(data),np.max(data)])
    width=interval[1]-interval[0]
    for i in range(int((1.-fraction)*len(data))): 
        integ=np.float64(0)
        j=np.argmin(np.abs(cdf-cdf[i]-np.float64(fraction)))
        if ((j < len(data)) & (j > i)):
            newwidth=data[j]-data[i]
            if newwidth < width:
                width=newwidth
                interval=[data[i],data[j]]

    return interval

def progressbar(percent):
    sys.stdout.write("%3d%%\r" % percent)
    sys.stdout.flush()

def readcdf(infile,niter):
    data=np.zeros(niter, np.float64)
    f=open(infile,'r')
    for i in range(niter):
        data[i]=np.float64(f.readline())
    datap,index=np.unique(data,return_index=True) #prune out repeated values
    cdf=index[1:]
    cdf=np.r_[cdf,len(data)]
    cdf=np.float64(cdf)/np.float64(niter) #convert counter to cdf
    return data,datap, cdf

def neatstats(infilepattern,niter,smooth=0.0075,avg=0,plots=0,verb=0):
    scale=np.linspace(0,1,num=101)
    ones=np.ones(101,dtype=np.float64)
    infilelist=glob.glob(infilepattern+'*_binned') 
    infilelist.sort() #see if you can find a better way of ordering them
    
    outfile=infilepattern+'_stats' #output file
    output=open(outfile,'w')
    #write stats header
    output.write("Quantity \t\t\t\t Average \t Confidence interval \t\t flags \n \n")

    i=-1
    for infile in infilelist:
        i+=1
        progressbar(100.*np.float(i)/np.float(len(infilelist)))
        infile=infile[:-7] #remove '_binned' from file name to get unbinned data

        #read data
        data, data_pruned, cdf =readcdf(infile,niter)

        #fit data with cubic spline
        fit=statsfit(data,data_pruned,cdf,niter,avg,smooth)

        #print
        outline='{:<45} {:>10.5g} {:>+10.5g} {:>+10.5g} {:>} \n'.format(infile,fit.average,fit.CI[1]-fit.average,fit.CI[0]-fit.average,fit.flags)
        output.write(outline)

        #output data if requested :)
        if verb > 0:
            pdffile=open(infile+'_cdfpdf','w')
            outline='{:>10.5g} {:>+10.5g} {:>+10.5g} {:>} \n'.format(fit.average,fit.CI[1]-fit.average,fit.CI[0]-fit.average,fit.flags)
            pdffile.write(outline)
            for j in range(len(data_prune)):
                outline='{:>10.5g} {:>+10.5g} {:>+10.5g}\n'.format(data_prune[j],cdf[i],fit.pdf[i])
                pdffile.write(outline)
            pdffline.close

        #begin plotting :)
        if plots > 0:
            if (plots == 2) | (len(fit.flags) > 0):
                plotfile=PdfPages(infile+'_plot.pdf')
                fig=plt.figure()
                cdfplot=fig.add_subplot(1,1,1)
                cdfplot.plot(data_pruned,cdf,color='blue')
                cdfplot.fill_between(data_pruned, cdf, 0, where=((data_pruned > fit.CI[0])&(data_pruned < fit.CI[1])), facecolor='blue', label='1 sigma region',alpha=0.35)
                cdfplot.plot(data_pruned[fit.iavg]*ones,cdf[fit.iavg]*scale,'--',color='black')
                cdfplot.plot(data_pruned[fit.iavg],cdf[fit.iavg],'o',color='black')
                plotfile.savefig()
                fig=plt.figure()
                pdfplot=fig.add_subplot(1,1,1)
                pdfplot.plot(data_pruned,fit.pdf,color='blue')
                pdfplot.fill_between(data_pruned, fit.pdf, 0, where=((data_pruned > fit.CI[0])&(data_pruned < fit.CI[1])), facecolor='blue', label='1 sigma region',alpha=0.35)
                pdfplot.plot(data_pruned[fit.iavg]*ones,fit.pdf[fit.iavg]*scale,'--',color='black')
                pdfplot.plot(data_pruned[fit.iavg],fit.pdf[fit.iavg],'o',color='black')
                plotfile.savefig()
                plotfile.close()
                plt.close('all')

    output.close
    print "Done."
    return 1 #replace with a proper status variable at some point :)

if __name__=="__main__":
    #use argparse to add lots of arguments and options :)
    parser=argparse.ArgumentParser(description='Reads the output from NEAT (using maxmimum verbosity) and attempts to approximate the cumulative distribution and probability density functions of the output quantities using spline fits. It then extracts an average (default: mode) and the smallest 68.3% confidence interval.')
    parser.add_argument('infilepattern',help='The name of the linelist analysed with NEAT. This can contain wildcards.',metavar='Linelist')
    parser.add_argument('niter',type=int,help='The number of Monte Carlo iterations used when running NEAT on the above linelist')
    parser.add_argument('-s','--smooth',default=0.0075,type=float,help='The value of the smoothing co-efficient to use in the spline fitting. Increasing this reduces oscillations in the derivatives of the spline, giving a smoother pdf. (Default: %(default)s)')
    parser.add_argument('-a','--average',default=0,type=int,choices=[0,1,2],help='Selects whether to use the mode (0), the median (1), or the mean (2) as the represenative average. If the mode has been chosen, but the pdf (and hence the mode) is not well defined it will default to the median. (Default: %(default)s)')
    parser.add_argument('-p','--plots',default=0,type=int,choices=[0,1,2],help='Controls whether plots of the CDF and PDF are generated. 0=No plots, 1=Only flagged distributions are plotted, 2=Plot everything. (Default: %(default)s)')
    parser.add_argument('-v','--verb',default=0,type=int,choices=[0,1],help='Controls output verbosity. For a value of 0, only a summary file is generated, while higher values produce ascii output of the cdf and fitted pdf in case you wish to generate your own plots of the distributions. (Default: %(default)s)')
    args=parser.parse_args()

    #arguments
    plots=args.plots
    smooth=args.smooth
    verb=args.verb
    avg=args.average
    infilepattern=args.infilepattern
    niter=args.niter

    a=neatstats(infilepattern,niter,smooth,avg,plots,verb)
