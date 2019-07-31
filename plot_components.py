import numpy as np
import matplotlib.pyplot as plt
from pylab import *
#import pyspeckit as ps
from scipy import io
from scipy import stats
from scipy.optimize import leastsq
#from lmfit import minimize, Parameters, Parameter, report_fit
#from lmfit.models import GaussianModel
import scipy.optimize as optimization
import matplotlib.ticker as ticker
import cmath as math
import pickle
import iminuit
import astropy.io.fits as pf
#from astropy.nddata import Cutout2D
import os,glob
#import string,math,sys,fileinput,glob,time
#load modules
#from pylab import *
import subprocess as sub

def get_ellipse_coords(a=0.0, b=0.0, x=0.0, y=0.0, angle=0.0, k=2):
    """ Draws an ellipse using (360*k + 1) discrete points; based on pseudo code
    given at http://en.wikipedia.org/wiki/Ellipse
    k = 1 means 361 points (degree by degree)
    a = major axis distance,
    b = minor axis distance,
    x = offset along the x-axis
    y = offset along the y-axis
    angle = clockwise rotation [in degrees] of the ellipse;
        * angle=0  : the ellipse is aligned with the positive x-axis
        * angle=30 : rotated 30 degrees clockwise from positive x-axis
    """
    pts = np.zeros((360*k+1, 2))

    beta = -angle * np.pi/180.0
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)
    alpha = np.radians(np.r_[0.:360.:1j*(360*k+1)])
 
    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    
    pts[:, 0] = x + (a * cos_alpha * cos_beta - b * sin_alpha * sin_beta)
    pts[:, 1] = y + (a * cos_alpha * sin_beta + b * sin_alpha * cos_beta)

    return pts

def ellipse_axis(x, y,s):
    x1=x-s
    x2=x+s
   
    if x1<x2:
        xaxis=np.linspace(x1,x2,50)    
    else:
        xaxis=np.linspace(x2,x1,50)           
    
    y1=y-s
    y2=y+s

    if y1<y2:
        yaxis=np.linspace(y1,y2,50)    
    else:
        yaxis=np.linspace(y2,y1,50)  

    return xaxis,yaxis


