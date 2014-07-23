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
    mean=np.zeros(1,dtype=np.float64)
#    imean=np.zeros(1,dtype=np.int)
    sigmamean=np.zeros(2,dtype=np.float64)
    sigmasymmetric=np.zeros(2,dtype=np.float64)
    sigmamedian=np.zeros(2,dtype=np.float64)
    nmodes=np.zeros(1,dtype=np.int)
    modes=np.zeros(1,dtype=np.float64)
    pdf=np.ndarray([])
    flags=""

def statsfit(data,cdf,niter,average,smooth=0.0075):
    #create object to contain correct values
    fits=fitdata()
    #first compute cubic spline fit
    splinefit=interpol.splrep(data,cdf,s=smooth)
    fits.CI=getci(data,cdf,0.683)
    #now find maxima
    firstderivative=interpol.splev(data,splinefit,der=1) 
    secondderivative=interpol.splev(data,splinefit,der=2)
    mode=data[np.argmax(firstderivative)]
    imode=np.argmax(firstderivative)
    fits.mode=mode
    fits.imode=imode
    fits.pdf=firstderivative
    fits.mean=np.mean(data)
    #check which roots are maxima & count them?
    
    #check which root has maximum value = mode

    #find 1-sigma interval
#    sigmamode=getinterval(data,cdf,mode,0.683)
#    fits.sigmamode=sigmamode
#    fits.sigmamean=getinterval(data,cdf,fits.mean,0.683)
    imed=niter/2 #check how python does int divison
    if (niter % 2) == 0: #find median
        median=(data[imed]+data[imed+1])/2.
    else:
        median=data[imed]
#    fits.sigmamedian=getinterval(data,cdf,median,0.683)
    fits.median=median
#    fits.sigmasymmetric=[data[imed - niter*0.3415],data[imed + niter*0.3415]]#getintervalsym(data,cdf,0.683)
    #if fits.sigmamode[0] == fits.mode:
    if fits.CI[0] == np.min(data):
        fits.flags+='l'
    #if fits.sigmamode[1] == fits.mode:
    if fits.CI[1] == np.min(data):
        fits.flags+='u'
    fits.average=fits.mode
    if average == 1:
        fits.average=fits.median
    if average == 2:
        fits.average=fits.mean
    if (average == 0) & (np.isnan(fits.pdf[0])):
        fits.average=fits.median
    #return only the important stuff :)
    return fits

def getintervalsym(data,cdf,imed,fraction):
    #integ=np.float64(0)
    #i=int(0)
    #while (integ < fraction) & (i < imed):
    #    integ=
    pass

def getinterval(data,cdf,point,fraction):
    #finds the smallest interval of data that contains fraction fraction of the cdf, and point point is inside it.
    interval=np.array([np.min(data),np.max(data)])
    #print interval
    width=interval[1]-interval[0]#len(data)
    ipoint=np.argmin(np.abs(data-point)) #get location of value that must lie within interval interval
    #print ipoint
    for i in range(ipoint+1): #range(len(data))
        integ=np.float64(0)
        #j=ipoint
        #while ((integ < fraction) & (j < len(data))):
        #   integ=cdf[j]-cdf[i]
        #   j=j+1
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
    #print interval
    width=interval[1]-interval[0]#len(data)
    #ipoint=np.argmin(np.abs(data-point)) #get location of value that must lie within interval interval
    #print ipoint
    for i in range(int((1.-fraction)*len(data))): #range(len(data))
        integ=np.float64(0)
        #j=ipoint
        #while ((integ < fraction) & (j < len(data))):
        #   integ=cdf[j]-cdf[i]
        #   j=j+1
        j=np.argmin(np.abs(cdf-cdf[i]-np.float64(fraction)))
        if ((j < len(data)) & (j > i)):
            newwidth=data[j]-data[i]#j-i#
            if newwidth < width:
                width=newwidth
                interval=[data[i],data[j]]

    return interval

def readcdf(infile,niter):
    data=np.zeros(niter, np.float64)
    cdf=np.float64(np.linspace(1./np.float64(niter),1.,niter,endpoint=True))
    f=open(infile,'r')
    for i in range(niter):
        data[i]=np.float64(f.readline())
    return data, cdf

if __name__=="__main__":
    #lookup argparse to add lots of arguments and options :)
    #including outline of argument parsing now so i don't forget how it works, even though it won't do anything useful just yet.
    parser=argparse.ArgumentParser(description='Reads the output from NEAT (using maxmimum verbosity) and attempts to approximate the cumulative distribution and probability density functions of the output quantities using spline fits. It then extracts an average (default: mode) and the smallest 68.3% confidence interval.')
    parser.add_argument('infilepattern',help='The name of the linelist analysed with NEAT. This can contain wildcards.',metavar='Linelist')
    parser.add_argument('niter',type=int,help='The number of Monte Carlo iterations used when running NEAT on the above linelist')
    parser.add_argument('-s','--smooth',default=0.0075,type=float,help='The value of the smoothing co-efficient to use in the spline fitting. Increasing this reduces oscillations in the derivatives of the spline, giving a smoother pdf. (Default: %(default)s)')
    parser.add_argument('-a','--average',default=0,type=int,choices=[0,1,2],help='Selects whether to use the mode (0), the median (1), or the mean (2) as the represenative average. If the mode has been chosen, but the pdf (and hence the mode) is not well defined it will default to the median. (Default: %(default)s)')
    parser.add_argument('-p','--plots',default=0,type=int,choices=[0,1,2],help='Controls whether plots of the CDF and PDF are generated. 0=No plots, 1=Only flagged distributions are plotted, 2=Plot everything. (Default: %(default)s)')
    parser.add_argument('-v','--verb',default=0,type=int,choices=[0,1,2],help='To be implemented. (Default: %(default)s)')
    args=parser.parse_args()

    plots=args.plots
    smooth=args.smooth
    verb=args.verb
    avg=args.average
    infilepattern=args.infilepattern
    niter=args.niter


    scale=np.linspace(0,1,num=101)
    ones=np.ones(101,dtype=np.float64)
    infilelist=glob.glob(infilepattern+'*_binned') #glob seems to create a strange order, but sorting isn't necessarily the most useful order, have to see if this can be improved
    infilelist.sort()
    
    outfile=infilepattern+'_stats' #output file
    output=open(outfile,'w')
    #write stats header
    output.write("Quantity \t\t\t\t Average \t Confidence interval \t\t flags \n \n")
#    texfike=outfile+'.tex'
    i=-1
    for infile in infilelist:
        i+=1
        infile=infile[:-7] #remove '_binned' from file name to get unbinned data
        print infile
        #read data
        data, cdf =readcdf(infile,niter)
        print np.max(data),np.min(data)
        #fit data with cubic spline
        fit=statsfit(data,cdf,niter,avg,smooth)
        #print
        outline='{:<45} {:>10.5g} {:>+10.5g} {:>+10.5g} {:>} \n'.format(infile,fit.average,fit.CI[1]-fit.average,fit.CI[0]-fit.average,fit.flags)
        output.write(outline)

        #begin plotting :)
        if plots > 0:
            if (plots == 2) | (len(fit.flags) > 0):
                print np.max(fit.pdf),np.min(fit.pdf)
                plotfile=PdfPages(infile+'_plot.pdf')
                fig=plt.figure()
                cdfplot=fig.add_subplot(1,1,1)
                cdfplot.plot(data,cdf,color='blue')
                cdfplot.fill_between(data, cdf, 0, where=((data > fit.sigmamode[0])&(data < fit.sigmamode[1])), facecolor='blue', label='1 sigma region',alpha=0.35)
                cdfplot.plot(data[fit.imode]*ones,cdf[fit.imode]*scale,'--',color='black')
                cdfplot.plot(data[fit.imode],cdf[fit.imode],'o',color='black')
                plotfile.savefig()
                fig=plt.figure()
                pdfplot=fig.add_subplot(1,1,1)
                pdfplot.plot(data,fit.pdf,color='blue')
                pdfplot.fill_between(data, fit.pdf, 0, where=((data > fit.sigmamode[0])&(data < fit.sigmamode[1])), facecolor='blue', label='1 sigma region',alpha=0.35)
                pdfplot.plot(data[fit.imode]*ones,fit.pdf[fit.imode]*scale,'--',color='black')
                pdfplot.plot(data[fit.imode],fit.pdf[fit.imode],'o',color='black')
                plotfile.savefig()
                plotfile.close()
                plt.close('all')

    output.close
    
