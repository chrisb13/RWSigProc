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
    plt.close('all')
    row=1
    col=2

    
    fig, axis = plt.subplots(\
        nrows=row, ncols=col, sharex=False, sharey=False,figsize=(5.5*col,10),\
        gridspec_kw={'hspace':.225,'wspace':.065}\
              )
    
    ax=axis[0]
    ax.grid(True)
    ax.contourf(df.columns,df.index,df.values,levels=levs,cmap=ccmap,extend='both')
    rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    rwsp.pl_inset_title_box(ax,'raw',bwidth="20%",location=1)

    ax=axis[1]
    ax.grid(True)
    size=11
    f=sc.array(sc.ones([size,size])/size**2)
    zf=nd.convolve(df.values,f)
    cs1=ax.contourf(df.columns,df.index,zf,levels=levs,cmap=ccmap,extend='both')
    rwsp.pl_inset_title_box(ax,'box filter',bwidth="40%",location=1)
    rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    plt.setp(ax.get_yticklabels(),visible=False)

    fig.colorbar(cs1, ax=axis.ravel().tolist(), pad=0.02, aspect = 30)

    fig.savefig(pout+'example_sqfilt_aviso.png',dpi=300,bbox_inches='tight')
    lg.info("Plot created: "+pout+'example_sqfilt_aviso.png')
    # plt.show()

    ###########################################################################
    #                          Do a sombrero filter                           #
    ###########################################################################

    plt.close('all')
    row=1
    col=2

    fig, axis = plt.subplots(\
        nrows=row, ncols=col, sharex=False, sharey=False,figsize=(5.5*col,10),\
        gridspec_kw={'hspace':.225,'wspace':.065}\
              )
    
    ax=axis[0]
    ax.grid(True)
    ax.contourf(df.columns,df.index,df.values,levels=levs,cmap=ccmap,extend='both')
    rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    rwsp.pl_inset_title_box(ax,'raw',bwidth="20%",location=1)

    ax=axis[1]
    ax.grid(True)
    f=rwsp.sombrero(11,11,-1000,365.25,sigma=2)
    zf=nd.convolve(df.values,f)
    cs1=ax.contourf(df.columns,df.index,zf,levels=levs,cmap=ccmap,extend='both')
    rwsp.pl_inset_title_box(ax,'sombrero filter',bwidth="40%",location=1)
    rwsp.change_tick_labels_add_dirs(ax,xyonly='x')
    plt.setp(ax.get_yticklabels(),visible=False)

    fig.colorbar(cs1, ax=axis.ravel().tolist(), pad=0.02, aspect = 30)

    fig.savefig(pout+'example_sombrero_aviso.png',dpi=300,bbox_inches='tight')
    lg.info("Plot created: "+pout+'example_sombrero_aviso.png')

    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
