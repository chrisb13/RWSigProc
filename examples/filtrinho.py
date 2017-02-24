#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 15:16:27 2017

@author: polito
"""

import numpy as np 
import scipy as sc
from scipy import ndimage as nd
import pylab as pl
z = np.random.normal(0,300,[100,300])
x, y = pl.meshgrid(pl.arange(1.,301,1)*7,pl.arange(1.,101,1)/4)
pl.subplot(211)
pl.contourf(x,y,z,pl.arange(-200.,201.,10.))
#pl.colorbar()

#eleven by eleven 'wide' and divided by 11^2 to ensure integrates to one (breaks the first rule)
f=sc.array(sc.ones([11,11])/11**2)

zf=nd.convolve(z,f)
pl.subplot(212)
pl.contourf(x,y,zf,pl.arange(-200.,201.,10.))
#pl.colorbar()

pl.show()



###########################################################################
#                             Proper attempt                              #
###########################################################################

#import recode as fir

#fir.sombrero()

##choose one wavelength by one period (in pixels/matrix elements)
#f=fir.sombrero(51,41,-1000,365.25,sigma=2)
#zf=nd.convolve(z,f)
