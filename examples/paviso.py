#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation: Climate Change Research Centre and ARC Centre of Excellence for Climate System Science.
#                Level 4, Mathews Building
#                University of New South Wales
#                Sydney, NSW, Australia, 2052
#   Contact: z3457920@unsw.edu.au
#   www:     christopherbull.com.au
#   Date created: Fri, 24 Feb 2017 18:37:12
#   Machine created on: ccrc165
#

"""
Example to filter and plot ./dt_global_allsat_madt_h_slice_single.nc 
"""
import sys,os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import collections
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from scipy import ndimage as nd
# yo! I changed this -----------
from scipy import signal as sig
# ------------------------------
import scipy as sc
from cb2logger import *

curdir=os.path.dirname(os.path.realpath(__file__))+'/'
print curdir

#insert path so we can import the package
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))+'/')

import RWSigProc as rwsp


if __name__ == "__main__": 
    LogStart('',fout=False)
    pout=curdir+'example_aviso_plots/'
    lg.info("Example of using AVISO, plots will be in: "+pout)

    rwsp.mkdir(curdir+'example_aviso_plots/')

    infile=curdir+'dt_global_allsat_madt_h_slice_single.nc'
    assert(os.path.exists(infile)),"AVISO file not found!"
    ifile=xr.open_dataset(infile)
    
    #munge data for pandas..
    data=collections.OrderedDict()
    for idx,ll in enumerate(ifile['lon'][:].values):
        data[ll]=ifile['adt'][:,0,idx]

    df=\
    pd.DataFrame(data,index=[pd.to_datetime(d) for d in ifile['time'][:].values] )

    df=df-df.mean()

    ###########################################################################
    #                      Do a 2D rectangular filter...                      #
    ###########################################################################

    ccmap='seismic'
    levs=np.linspace(-0.3,0.3,30)
    # plt.close('all')
    # row=1
    # col=2

    
    # fig, axis = plt.subplots(\
    #     nrows=row, ncols=col, sharex=False, sharey=False,figsize=(5.5*col,10),\
    #     gridspec_kw={'hspace':.225,'wspace':.065}\
    #           )
    
    # ax=axis[0]
    # ax.grid(True)
    # ax.contourf(df.columns,df.index,df.values,levels=levs,cmap=ccmap,extend='both')
    # rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    # rwsp.pl_inset_title_box(ax,'raw',bwidth="20%",location=1)

    # ax=axis[1]
    # ax.grid(True)
    # size=11
    # f=sc.array(sc.ones([size,size])/size**2)
    # zf=nd.convolve(df.values,f)
    # cs1=ax.contourf(df.columns,df.index,zf,levels=levs,cmap=ccmap,extend='both')
    # rwsp.pl_inset_title_box(ax,'box filter',bwidth="40%",location=1)
    # rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    # plt.setp(ax.get_yticklabels(),visible=False)

    # fig.colorbar(cs1, ax=axis.ravel().tolist(), pad=0.02, aspect = 30)

    # fig.savefig(pout+'example_sqfilt_aviso.png',dpi=300,bbox_inches='tight')
    # lg.info("Plot created: "+pout+'example_sqfilt_aviso.png')
    # plt.show()

    ###########################################################################
    #                          Do a sombrero filter                           #
    ###########################################################################

    plt.close('all')
    row=1
    col=5

    fig, axis = plt.subplots(\
        nrows=row, ncols=col, sharex=False, sharey=False,figsize=(5.5*col,10),\
        gridspec_kw={'hspace':.225,'wspace':.065}\
              )
    
    ax=axis[0]
    ax.grid(True)
    ax.contourf(df.columns,df.index,df.values,levels=levs,cmap=ccmap,extend='both')
    rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    rwsp.pl_inset_title_box(ax,'raw',bwidth="20%",location=1)

    # this block is for the annual signal, non-propagating, lambda gigantic

    fa = rwsp.sombrero(365,201,-1.e9,365.25,sigma=2)
#    zf=nd.convolve(df.values,fa) # memory issues with large matrices, try this:
    za = sig.fftconvolve(df.values,fa,mode='same')
    za = za*rwsp.maxvar(df.values,za)
    # after you are happy with the filter result, remove it and look for
    # the next dominant signal out of the remainder (zrem)
    zrem = df.values-za

    # this block is for the annual Rossby waves,

    fr1 = rwsp.sombrero(365,61,-60.,365.25,sigma=2)
    zr1 = sig.fftconvolve(zrem,fr1,mode='same')
    zr1 = zr1*rwsp.maxvar(zrem,zr1)

    zrem = zrem-zr1

    # this block is for semiannual Rossby waves

    fr2 = rwsp.sombrero(183,31,-30.,183.,sigma=2)
    zr2 = sig.fftconvolve(zrem,fr2,mode='same')
    zr2 = zr2*rwsp.maxvar(zrem,zr2)

    zrem = zrem-zr2
 

    ax=axis[1]
    ax.grid(True)
    cs1=ax.contourf(df.columns,df.index,za,levels=levs,cmap=ccmap,extend='both')
    rwsp.pl_inset_title_box(ax,'Annual signal',bwidth="40%",location=1)
    rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    plt.setp(ax.get_yticklabels(),visible=False)

    ax=axis[2]
    ax.grid(True)
    cs2=ax.contourf(df.columns,df.index,zr1,levels=levs,cmap=ccmap,extend='both')
    rwsp.pl_inset_title_box(ax,'12mo Rossby',bwidth="40%",location=1)
    rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    plt.setp(ax.get_yticklabels(),visible=False)

    ax=axis[3]
    ax.grid(True)
    cs3=ax.contourf(df.columns,df.index,zr2,levels=levs,cmap=ccmap,extend='both')
    rwsp.pl_inset_title_box(ax,'6mo Rossby',bwidth="40%",location=1)
    rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    plt.setp(ax.get_yticklabels(),visible=False)

    ax=axis[4]
    ax.grid(True)
    cs4=ax.contourf(df.columns,df.index,zrem,levels=levs,cmap=ccmap,extend='both')
    rwsp.pl_inset_title_box(ax,'Unexplained',bwidth="40%",location=1)
    rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    plt.setp(ax.get_yticklabels(),visible=False)

    fig.colorbar(cs1, ax=axis.ravel().tolist(), pad=0.02, aspect = 30)

    fig.savefig(pout+'example_sombrero_aviso.png',dpi=300,bbox_inches='tight')
    lg.info("Plot created: "+pout+'example_sombrero_aviso.png')

    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
