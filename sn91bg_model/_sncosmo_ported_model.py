#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module defines the SN91bgSource class, which acts as a custom SNCosmo
source object for 91bg like supernovae.
"""

import os

import numpy as np
import sncosmo
from scipy.interpolate import RectBivariateSpline,interpn
import matplotlib.pyplot as plt

FILE_DIR = os.path.abspath(os.path.dirname(__file__))
COMPILED_MODEL_PATH = os.path.join(FILE_DIR, 'complete_template.npy')

def bi_search(a,x):
    """
    Binary search

    Args:
        a   (list): the list in which the number x will be searched
        x    (num): the number to be searched 

    Returns:
        the position of nearest neighbors of the number x
    """
    if x in a:
        return a.index(x)
    l,r=0,len(a)
    while(abs(r-l)>1):
        m=(l+r)//2
        if a[m]>x:
            r=m
        if a[m]<x:
            l=m
    return l,r

def bilinear_interp(x1,x2,y1,y2,f,x,y):
    """
    Bilinear interpolate point (x,y) in the grid (x1,y1), (x1,y2), (x2,y1), (x2,y2)

    Args:
        x1,x2,y1,y2  (num): grid point coordinates
        x,y          (num): coordinate of the interpolating point
        f           (list): a list containing four values at (x1,y1), (x1,y2), (x2,y1), (x2,y2)

    Returns:
        the interpolating value
    """
    f11=f[0]
    f12=f[1]
    f21=f[2]
    f22=f[3]
    fy1=(x2-x)*f11/0.1+(x-x1)*f21/0.1
    fy2=(x2-x)*f12/0.1+(x-x1)*f22/0.1
    f_interp=(y2-y)*fy1/0.25+(y-y1)*fy2/0.25
    return f_interp

class SN91bgSource(sncosmo.Source):
    _param_names = ['amplitude', 'stretch', 'color']
    param_names_latex = ['A', 'x1', 'c']   # used in plotting display

    def __init__(self):
        super(SN91bgSource, self).__init__()

        self.name = '91bg model'
        self.version = 'None'

        # Model spectra and initial guess for stretch, color, and amplitude
        self._flux_values = self._get_91bg_model()
        self._parameters = np.array([1., 1., 0.55])

        # 4-dimension grid points in the model
        self._stretch = np.array([0.65,0.75,0.85,0.95,1.05,1.15,1.25])
        self._color = np.array([0.0, 0.25, 0.5, 0.75, 1])
        self._phase = np.arange(-13, 101, 1)
        self._wave = np.arange(1000, 12001, 10)

        # Creat bi-cubic spline in dimension phase and wavelength for all 35 templates
        self._model_flux=np.full((len(self._stretch),len(self._color)),None).tolist()
        for i in range(len(self._stretch)):
            for j in range(len(self._color)):
                self._model_flux[i][j]=RectBivariateSpline(self._phase, self._wave, self._flux_values[i][j], kx=3, ky=3)

    @staticmethod
    def _get_91bg_model():
        """Load template spectra for SN 1991bg

        Template spans 7 stretch values, 5 color values, 114 dates,
        and 1101 wavelengths.

        Returns:
            An array of shape (7, 5, 114, 1101)
        """

        return np.load(COMPILED_MODEL_PATH)

    def _flux(self, phase, wave):
        """Bilinear Interpolate template flux for stretch, color

        Time and wavelength are given by arguments phase and wave. 
        Stretch, color, and amplitude are given by self._parameters

        Args:
            phase (list): A list of days till maximum for the desired flux
            wave  (list): A list of wavelengths for the desired flux

        Returns:
            A 2d array of flux with shape (<len(phase)>, <len(wave)>)
        """

        A, x1, c=self._parameters

        # Find the 4 nearest grid points of (x1,c)
        st1,st2=bi_search(self._stretch,x1)
        c1,c2=bi_search(self._color,c)
        f=[self._model_flux[st1][c1](phase,wave),self._model_flux[st1][c2](phase,wave),\
          self._model_flux[st2][c1](phase,wave),self._model_flux[st2][c2](phase,wave)]

        # Bilinear interpolation for stretch and color
        return A*bilinear_interp(self._stretch[st1],self._stretch[st2],self._color[c1],self._color[c2],f,x1,c)