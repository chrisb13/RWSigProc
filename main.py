#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation: Climate Change Research Centre and ARC Centre of Excellence for Climate System Science.
#                Level 4, Mathews Building
#                University of New South Wales
#                Sydney, NSW, Australia, 2052
#   Contact: z3457920@unsw.edu.au
#   www:     christopherbull.com.au
#   Date created: Fri, 24 Feb 2017 16:54:19
#   Machine created on: chris-VirtualBox2
#
"""
Filter module file
"""

import numpy as np
import matplotlib.pyplot as plt
        

def sombrero(v,h,L,T,sigma=2):
    """produces a sombrero-like shaped matrix z, used for FIR filtering
    
    :v: rows
    :h: columns
    :L: L horizontal wavelength in matrix units, may be negative, default v
    :T: T vertical wavelength in matrix units, positive, default h

    Notes
    -------
    Try and make v and h odd numbers!
    
    :returns: array that is v rows by h columns
    """
    #% get the gaussian envelope

    s=sigma  #% s controls the value at the borders
         #% s is the sigma of the normal distribution

    #[x,y]=meshgrid(0.5:h-.5,0.5:v-.5);
    #x=(pi*(x-h/2)/(h*s)).^2; y=(pi*(y-v/2)/(v*s)).^2;
    #g=exp(-0.5*(x+y)) ./ (sqrt(2*pi) .* s);
    #g=g-min(g(:));
    #g=g/max(g(:));

    x,y=np.meshgrid(np.arange(0.5,h+0.5,1),np.arange(0.5,v+0.5,1))
    #print '--'
    #print x
    #print '--'
    #print y
    #print '--'
    
    #x=(pi*(x-h/2)/(h*s)).^2; y=(pi*(y-v/2)/(v*s)).^2;

    x=(np.pi*(x-h/2)/(h*s))**2
    y=(np.pi*(y-v/2)/(v*s))**2

    g=np.exp(-0.5*(x+y)) / (np.sqrt(2*np.pi) * s)

    g=g-g.min()
    g=g/g.max()

    #% get the cosine surface and multiply by the gaussian

    #[x,y]=meshgrid( ((1:h)-(h+1)/2) , ((1:v)-(v+1)/2) );
    x,y=np.meshgrid((np.arange(h)+1)-(h+1)/2,(np.arange(v)+1)-(v+1)/2)

    #print '--'
    #print x
    #print '--'
    #print y
    #print '--'
    
    z=g*np.cos((2*np.pi/L)*x - (2*np.pi/T)*y)

    #% accept negative values, normalize for zero sum
    sz=np.sum(z[:])
    s1=np.sum(np.ones(np.size(z)))
    z=z-sz/s1
    sa=np.sum(np.abs(z[:]))
    z=z/sa
    return z

if __name__ == "__main__": 
    #LogStart('',fout=False)
    from _cblogger import _LogStart
    _lg=_LogStart().setup()
    import time

    #choose one wavelength by one period (in pixels/matrix elements)
    #annual #(Big wavelength!)
    z=sombrero(51,41,-1000,365.25,sigma=2)
    plt.contourf(z,30)
    plt.show()

    _lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    _lg.info("Local current time : "+ str(localtime))
    _lg.info('SCRIPT ended')

else:
    from ._cblogger import _LogStart
    _lg=_LogStart().setup()
