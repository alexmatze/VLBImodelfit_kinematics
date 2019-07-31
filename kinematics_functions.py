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
import os,glob
#import string,math,sys,fileinput,glob,time
#load modules
#from pylab import *
import subprocess as sub
from plot_components import get_ellipse_coords, ellipse_axis
import urllib2
from astropy import units as u
#from astropy.coordinates import SkyCoord

#FUNCTION TO READ THE HEADER AND TAKE IMPORTANT PARAMETERS AS
#cell
#BMAJ, BMIN, BPA
#date, freq and epoch

def axis_limits(y):
	ymin = np.min(y)
	if -1 <= ymin <= 1: 	
		ymin_lim = ymin - 0.1
	else:
		ymin_lim = ymin - 0.2*ymin
	ymax = np.max(y)
	if -1 <= ymin <= 1: 	
		ymax_lim = ymax + 0.1
	else:
		ymax_lim = ymax + 0.2*ymax

	return ymin_lim,ymax_lim


def search_rms():
	search_string = 'Estimated noise='
	rms = []

	with open('difmap.log','r') as infile:
		for line in infile:
			if search_string in line:
				s=line
				s1 = s.split()
				s2 = s1[2].split('=')
				rms_value = float(s2[1])
				if s1[3][0]== 'm':
					rms.append(rms_value/1000.)
					
				elif s1[3][0] == 'J':
					rms.append(rms_value)

	return rms


def take_header(x):
	'''
	given x MAP file to be read
	'''

	#open the fits file and read the header
	FITS = pf.open(x)
	HEADER = pf.getheader(x)

	#read the parameters that contain the cell, bmaj, bmin and bpa 
	cell = np.abs(HEADER['CDELT1']*3600*1000)  #= -9.7222222253244E-08 /Pixel increment        
	BMAJ =HEADER['BMAJ']*3600*1000
	BMIN = HEADER['BMIN']*3600*1000  # /Reference pixel    
	BPA = HEADER['BPA']

	#read the date as it is in the header
	date = HEADER['DATE-OBS']

	#for dates in the format year-month-day
	if date[2].isdigit():
		d=date.split('-')
		date=str(d[0])+'_'+str(d[1])+'_'+str(d[2])
		epoch=(float(d[0])+((float(d[1])-1)*30.4375+float(d[2]))/365.25)
	#for dates in the format day/month/year (with only the two last digits of the year)
	else:
		d=date.split('/')
		print d[2][0]
		if d[2][0]=='9':
                	d[2]=1900+float(d[2])
                	date=str("%0.0f" % (d[2],))+'_'+str(d[1])+'_'+str(d[0])
                	epoch=(float(d[2])+((float(d[1])-1)*30.4375+float(d[0]))/365.25)        
		else:
                	d[2]=2000+float(d[2])        
                	date=str(d[2])+'_'+str(d[1])+'_'+str(d[1])
                	epoch=(float(d[2])+((float(d[1])-1)*30.4375+float(d[0]))/365.25)
 	freq = float(HEADER['CRVAL3'])/1e9

	#equivalent cicular beam
	circ_beam = np.sqrt(BMAJ*BMIN)

	source_name = HEADER['OBJECT']

	rms = HEADER['NOISE']

	mapsize =HEADER['NAXIS1'] #mapsize           
	datamax = HEADER['DATAMAX']   

	return cell, BMAJ,BMIN,BPA,date,freq,epoch,circ_beam, source_name,rms,mapsize,datamax
	#mapsize [10], pixelsize[0]

def order_by_date(files,models,modelfit,fits):

	#initialize arrays
	epoch = np.array([0.]*len(fits))
	cell = np.array([0.]*len(fits))
	bmaj = np.array([0.]*len(fits))
	bmin = np.array([0.]*len(fits))
	bpa = np.array([0.]*len(fits))
	beam = np.array([0.]*len(fits))
	size_map = np.array([0.]*len(fits))
	size_map_y =np.array([0.]*len(fits))
	
	#read some parameters of the map
	for i in xrange(0,len(fits)):
		header = take_header(fits[i])
		epoch[i] = header[6]
		cell[i] = header[0]
		bmaj[i] = header[1]
		bmin[i] = header[2]
		bpa[i] = header[3]
		beam[i] = header[7]
		
	#store them in ascendent order related with the date
	n=0
	for z in epoch:
		n+=1
	for i in range(n):
		for j in range(1,n-i):
			if epoch[j-1] > epoch[j]:
				(epoch[j-1],epoch[j]) = (epoch[j],epoch[j-1])
				(cell[j-1],cell[j]) = (cell[j],cell[j-1])
				(bmaj[j-1],bmaj[j]) = (bmaj[j],bmaj[j-1])
				(bmin[j-1],bmin[j]) = (bmin[j],bmin[j-1])
				(bpa[j-1],bpa[j]) = (bpa[j],bpa[j-1])
				(beam[j-1],beam[j]) = (beam[j],beam[j-1])
				(size_map[j-1],size_map[j]) = (size_map[j],size_map[j-1])
				(size_map_y[j-1],size_map_y[j]) = (size_map_y[j],size_map_y[j-1])
				(files[j-1],files[j]) = (files[j],files[j-1])
				(models[j-1],models[j]) = (models[j],models[j-1])
				(fits[j-1],fits[j]) = (fits[j],fits[j-1])
				(modelfit[j-1],modelfit[j]) = (modelfit[j],modelfit[j-1])

	return epoch,cell,bmaj,bmin,bpa,beam,files,models,fits,modelfit

def read_modfile(file1,beam,bmaj,bmin,z,freq,noise,errors):

	temp_filenoise=np.loadtxt('rmsvalue.dat')
	noise = temp_filenoise/1000.

	nfiles = len(file1)

	r_arr = []
	errr_arr = [] #np.array([0.]*nfiles)
	psi_arr =  []
	errpsi_arr =  []
	size_arr =  []
	errsize_arr =  []
	flux_arr =  []
	errflux_arr =  []
	tb_arr =  []
	errtb_arr =  []

	cte = 1.2*10**9

	ntot=0
	for k in xrange (0,nfiles):	
		with open(file1[k]) as myfile:
		    	count = sum(1 for line in myfile if line.rstrip('\n'))
	
		count = count-4

		#n = len(rms[k])

		n = count

		split_f=[]
		c=[]
		r=np.array([0.]*n)
		errr=np.array([0.]*n)
		psi=np.array([0.]*n)
		errpsi=np.array([0.]*n)
		size=np.array([0.]*n)
		errsize=np.array([0.]*n)
		tb=np.array([0.]*n)
		errtb=np.array([0.]*n)
		tbNew=np.array([0.]*n)
		errtbNew=np.array([0.]*n)
		flux=np.array([0.]*n)
		fluxpeak = np.array([0.]*n)
		rms = np.array([0.]*n)
		errflux=np.array([0.]*n)
		lim_resol=np.array([0.]*n)
		errlim_resol=np.array([0.]*n)
		
		temp=file1[k]
		temp_file=open(temp,mode='r')
		temp_file.readline()
		temp_file.readline()
		temp_file.readline()
		temp_file.readline()
		for i in xrange(0,n):
        		split_f = temp_file.readline().split()
        		flux[i] = (float(split_f[0][:-1]))
        		r[i] = (float(split_f[1][:-1]))
        		psi[i] = (float(split_f[2][:-1])*np.pi/180.)
        		size[i] = (float(split_f[3][:-1])/2.)
        		tb[i] = (float(split_f[7])) 

		if errors == True:
			temp_file2=open('pos_errors.dat',mode='r')
			temp_file2.readline()
			temp_file2.readline()
			for i in xrange(0,ntot):
				temp_file2.readline()
			for i in xrange(0,n):		
				split_f = temp_file2.readline().split()
				fluxpeak[i] = (float(split_f[2][:-1]))	
				rms[i] = (float(split_f[1][:-1]))			

			for i in xrange(0,n):
				errflux[i] = rms[i]
				snr = fluxpeak[i]/rms[i]#[k][i] #change to flux_peak

				dlim = 4/np.pi*np.sqrt(np.pi*np.log(2)*beam*np.log((snr)/(snr-1.))) #np.log((snr+1.)/(snr))) 4/np.pi*beam
				if size[i] > beam:
					ddec=np.sqrt(size[i]**2-beam**2)
				else:
					ddec=0.

				y=[dlim,ddec]
				dg=np.max(y)   
				err_size = rms[i]*dlim/fluxpeak[i]
				err_r = err_size/2.
				if r[i] > 0.:
					err_psi = np.real(math.atan(err_r*180./(np.pi*r[i])))
				else:
					err_psi = 1./5*beam
				
				if err_size < 2./5*beam:	
					errsize[i] = 2./5*beam
				else:
					errsize[i] = (err_size)
				if err_r < 1./5*beam:	
					errr[i] =  1./5*beam
					if errr[i] < 1/2*size[i]:
						errr[i] = 1/2*size[i]
				else:
					errr[i] = (err_r)
				errpsi[i] = (err_psi)	
				if dlim < beam:
					tbNew[i] = cte*(1.+z)/(freq**2)*fluxpeak[i]/(dlim**2)
					errtbNew[i] = tbNew[i]*np.sqrt((rms[i]/fluxpeak[i])**2)
				else:
					tbNew[i] = cte*(1.+z)/(freq**2)*fluxpeak[i]/(size[i]**2)
					errtbNew[i] = tbNew[i]*np.sqrt((rms[i]/fluxpeak[i])**2+(4*errsize[i]/size[i])**2)

		elif errors == 'Done':
			print 'done'

			#errr = size * (3 * rms)/flux
			#and if that is smaller than 1/5 the beam, put that as lower limit
			#rms and flux have to be in the same unit

		else:
			for i in xrange(0,n):
				errflux[i] = 0.1*flux[i]
				#print 'here'
                                #print 'size,bmaj,snr', size[i],bmaj[k],flux[i]/(3*noise[k]),flux[i],3*noise[k]
				#if size[i] < bmaj[k]:
				#	errr[i] = np.max([bmaj[k]*(3*noise[k])/flux[i],1./5*bmaj[k],1./2*size[i]])
				#else:
				#	errr[i] = np.max([size[i]*(3*noise[k])/flux[i]*size[i]**2/(bmaj[k]*bmin[k]),1./5*bmaj[k],1./2*size[i]])
				#errpsi[i] = 0.	
				#errsize[i] = 2./5*beam	
                                errsize[i] = 3*noise[k]*size[i]/flux[i]
                                errr[i] = errsize[i]/2.

				if r[i] > 0.:
					err_psi = np.real(math.atan(180.*errr[i]/(np.pi*r[i])))
				else:
					err_psi = 1./5*beam	
				
				print np.rad2deg(psi[i]),np.rad2deg(err_psi)
                                

				errpsi[i] = np.rad2deg(err_psi)		
                print np.mean(errpsi),np.std(errpsi)
		r_arr.append(r)
		errr_arr.append(errr)
		psi_arr.append(psi)
		errpsi_arr.append(errpsi)
		size_arr.append(size)
		errsize_arr.append(errsize)
		flux_arr.append(flux)
		errflux_arr.append(errflux)
		if errors == True:
			tb_arr.append(tbNew)
			errtb_arr.append(errtbNew)
		else:
			tb_arr.append(tb)
			errtb_arr.append(errtbNew)
	
		ntot = n + ntot + 1

	return r_arr,errr_arr,psi_arr,errpsi_arr,size_arr,errsize_arr,tb_arr,errtb_arr,flux_arr,errflux_arr

def x_y(r,errr,psi,errpsi,errors):
	n = len(r)
	x,errx = np.array([0.]*n),np.array([0.]*n)
	y,erry = np.array([0.]*n),np.array([0.]*n)
	x_arr, errx_arr = [], []
	y_arr, erry_arr = [], []
	for i in xrange (0,n):
		x=r[i]*np.sin(psi[i])
	        y=r[i]*np.cos(psi[i])
		if errors == True:
			errx=np.sqrt((errr[i]*np.cos(psi[i]))**2+(r[i]*np.sin(psi[i])*errpsi[i])**2)
			erry=np.sqrt((errr[i]*np.sin(psi[i]))**2+(r[i]*np.cos(psi[i])*errpsi[i])**2)
		else:
			errx = errr[i]	
			erry = errr[i]

		x_arr.append(x)
		errx_arr.append(errx)
		y_arr.append(y)
		erry_arr.append(erry)

	x_arr = np.asarray(x_arr)
	errx_arr = np.asarray(errx_arr)
	y_arr = np.asarray(y_arr)
	erry_arr = np.asarray(erry_arr)

    	return x_arr,errx_arr,y_arr,erry_arr

def r_psi(x,errx,y,erry):
	n = len(r)
	r,errr = np.array([0.]*n),np.array([0.]*n)
	psi,errpsi = np.array([0.]*n),np.array([0.]*n)
	r_arr, errr_arr = [], []
	psi_arr, errpsi_arr = [], []
	for i in xrange (0,n):
		r=np.sqrt(x[i]**2+y[i]**2)
	        psi=np.atan(y[i]/x[i])
		#errr=np.sqrt((1/(2*r)*2*x[i]*errx[i])**2+(1/(2*r)*2*y[i]*erry[i])**2)
		#errpsi=np.sqrt(((y[i]/([x[i]**2+y[i])**2])*errx[i])**2+((x[i]/([x[i]**2+y[i])**2])*erry[i])**2)

		r_arr.append(r)
		#errr_arr.append(errr)
		psi_arr.append(psi)
		#errpsi_arr.append(errpsi)

    	return r_arr,psi_arr

def find_nearest(array,value):
	index = (np.abs(array-value)).argmin()
	return array[index], index

def convolve_difmap(files,models,bmaj,bmin,bpa,xshift,yshift,mapsize,pixelsize,files2):
	difmap_location = sub.call(["which difmap"], shell=True)
	#open difmap and convolve maps with new beam
	for i in xrange(0,len(files)): #/usr/local/pgplot.old/uvf_difmap/difmap 
		p=sub.Popen(['/usr/local/bin/difmap'],stdin=sub.PIPE,stdout=sub.PIPE,close_fds=True,universal_newlines=True)
		(child_stdout,child_stdin)=(p.stdout,p.stdin)
		child_stdin.write('observe %s\n' %files[i])
		child_stdin.write('select PI\n')
		child_stdin.write('uvweight %s,%s\n'%(0,-1)) #WARNING read it from the par ####
		child_stdin.write('mapsize %s,%s\n' %(mapsize,pixelsize))
		child_stdin.write('rmod %s\n' %models[i])
		child_stdin.write('restore %s,%s,%s\n' %(bmaj,bmin,bpa))
		if float(xshift)!=0 and float(yshift)!=0:
			child_stdin.write('shift %3.3f, %3.3f\n' %(xshift,yshift))
		child_stdin.write('wmap %s\n' %files2[i])
		child_stdin.write('quit\n')
		p.wait()
	
	#remove log file
	os.system('rm difmap.log*\n')

def shift_modelfit_difmap(files,modelfit_mods,bmaj,bmin,bpa,xshift,yshift,mapsize,pixelsize,files2):
	#open difmap and convolve maps with new beam
	from distutils.spawn import find_executable
	difmap_path = find_executable('difmap')
	for i in xrange(0,len(files)): #/usr/local/pgplot.old/uvf_difmap/difmap 
		path = difmap_path #'/usr/local/bin/difmap'
		p=sub.Popen([path],stdin=sub.PIPE,stdout=sub.PIPE,close_fds=True,universal_newlines=True)
		(child_stdout,child_stdin)=(p.stdout,p.stdin)
		child_stdin.write('observe %s\n' %files[i])
		child_stdin.write('select PI\n')
		child_stdin.write('uvweight %s,%s\n'%(0,-1)) #WARNING read it from the par ####
		child_stdin.write('mapsize %s,%s\n' %(mapsize,pixelsize))
		child_stdin.write('rmod %s\n' %modelfit_mods[i])
		child_stdin.write('restore %s,%s,%s\n' %(bmaj,bmin,bpa))
		if float(xshift)!=0 and float(yshift)!=0:
			child_stdin.write('shift %3.3f, %3.3f\n' %(xshift,yshift))
		child_stdin.write('wmod %s\n' %files2[i])
		child_stdin.write('quit\n')
		p.wait()
	
	#remove log file
	os.system('rm difmap.log*\n')

def ellipse_axis_lines(x,y,size):
	pts_arr=[]
	pt_arr=[]
	x_el_arr=[]
	x_elH_arr=[]
	y_el_arr=[]
	y_elH_arr=[]

	for i in xrange(0,len(x)):
		n = len(x[i])
		pts, pt = [], []
		x_el, y_el = [], []
		x_elH, y_elH = [], []
		for k in xrange(0,n):
			pts.append(get_ellipse_coords(a=size[i][k], b=size[i][k], x=x[i][k],y=y[i][k], angle=0))
			pt.append(get_ellipse_coords(a=0.01, b=0.01, x=x[i][k],y=y[i][k], angle=0))
			#lines axis ellipses      
			x_el.append(ellipse_axis(x=float(x[i][k]),y=float(y[i][k]),s=float(size[i][k]))[0])
			y_el.append(ellipse_axis(x=x[i][k],y=y[i][k],s=size[i][k])[1])
			x_elH.append(np.linspace(x[i][k],x[i][k],50))
			y_elH.append(np.linspace(y[i][k],y[i][k],50))

		pts_arr.append(pts)
		pt_arr.append(pt)
		x_el_arr.append(x_el)
		y_el_arr.append(y_el)
		x_elH_arr.append(x_elH)
		y_elH_arr.append(y_elH)

	return pts_arr,pt_arr,x_el_arr,y_el_arr,x_elH_arr,y_elH_arr

def plot_components(pts_arr,x_el_arr,x_elH_arr,y_elH_arr,y_el_arr,arr,orientation):
    	for i in xrange(0,len(x_el_arr)):
		arr_n = arr[i]
		for j in xrange(0,len(x_el_arr[i])):
			if orientation == 'h':
				plt.plot(pts_arr[i][j][:,0], pts_arr[i][j][:,1]-arr_n, color='blue',linewidth=4)
				plt.plot(x_el_arr[i][j], y_elH_arr[i][j]-arr_n, color='blue',linewidth=4) 
				plt.plot(x_elH_arr[i][j], y_el_arr[i][j]-arr_n, color='blue',linewidth=4)
			if orientation == 'v':
				plt.plot(pts_arr[i][j][:,0]-arr_n, pts_arr[i][j][:,1], color='blue',linewidth=4)
				plt.plot(x_el_arr[i][j]-arr_n, y_elH_arr[i][j], color='blue',linewidth=4) 
				plt.plot(x_elH_arr[i][j]-arr_n, y_el_arr[i][j], color='blue',linewidth=4)

def plot_maps(realDAT,ext,arr,first_contour,orientation):
	for i in xrange(0,len(realDAT)):
		arr_n = arr[i]
		levels = first_contour[i]*realDAT[i].std()*np.array([-1., 1., 1.41,2.,2.83,4.,5.66,8.,11.3,16.,
	                                22.6,32.,45.3,64.,90.5,128.,181.,256.,362.,512.,
	                                724.,1020.,1450.,2050.])
		if orientation == 'h':
			cset = plt.contour(realDAT[i], levels, inline=1,
		                  colors=['grey']
		                  ,extent=[ext[0],ext[1],ext[2]-arr_n,ext[3]-arr_n], aspect=1.0
		                  )
		if orientation == 'v':
			cset = plt.contour(realDAT[i], levels, inline=1,
		                  colors=['grey']
		                  ,extent=[ext[0]-arr_n,ext[1]-arr_n,ext[2],ext[3]], aspect=1.0
		                  )

def plot_related_components(pts_arr,x_el_arr,x_elH_arr,y_elH_arr,y_el_arr,arr,index,color_comp,orientation):
	for i in xrange(0,len(x_el_arr)):
		arr_n = arr[i]
		idx = int(index[i])
		if orientation == 'h':
			plt.plot(pts_arr[i][idx][:,0], pts_arr[i][idx][:,1]-arr_n, color=color_comp,linewidth=4)
	       		plt.plot(x_el_arr[i][idx], y_elH_arr[i][idx]-arr_n, color=color_comp,linewidth=4) 
	       		plt.plot(x_elH_arr[i][idx], y_el_arr[i][idx]-arr_n, color=color_comp,linewidth=4)
		if orientation == 'v':
			plt.plot(pts_arr[i][idx][:,0]-arr_n, pts_arr[i][idx][:,1], color=color_comp,linewidth=4)
	       		plt.plot(x_el_arr[i][idx]-arr_n, y_elH_arr[i][idx], color=color_comp,linewidth=4) 
	       		plt.plot(x_elH_arr[i][idx]-arr_n, y_el_arr[i][idx], color=color_comp,linewidth=4)

def plot_NiA_components(pts,xel,yelH,xelH,yel,arr,indexComp,colorComp,lastComp,ptstemp,xeltemp,yelHtemp,xelHtemp,yeltemp,arrtemp,indexComptemp,orientation):
# plot_NiA_components(KinematicsWindow.WINpts_arrNiA,KinematicsWindow.WINx_el_arrNiA,KinematicsWindow.WINy_elH_arrNiA,
#			KinematicsWindow.WINx_elH_arrNiA,KinematicsWindow.WINy_el_arrNiA,KinematicsWindow.WINarrNiA,
#			KinematicsWindow.WINcompAutomaticIndexNotInAll,color_comp,last_comp,
#			KinematicsWindow.WINpts_arrNiAtemp,KinematicsWindow.WINx_el_arrNiAtemp,KinematicsWindow.WINy_elH_arrNiAtemp,
#			KinematicsWindow.WINx_elH_arrNiAtemp,KinematicsWindow.WINy_el_arrNiAtemp,KinematicsWindow.WINarrNiAtemp,
#			KinematicsWindow.WINcompAutomaticIndexNiAnewComp
	for i in xrange(0,len(xel)):
		arrn= arr[i]
		idxList = indexComp[i]
		for j in xrange(0,len(xel[i])):
			arr_n = arrn[j]
			idx = int(idxList[j])
			if orientation == 'h':
				plt.plot(pts[i][j][:,0], pts[i][j][:,1]-arr_n, color=colorComp[lastComp+i+1],linewidth=4)
				plt.plot(xel[i][j], yelH[i][j]-arr_n, color=colorComp[lastComp+i+1],linewidth=4) 
				plt.plot(xelH[i][j], yel[i][j]-arr_n, color=colorComp[lastComp+i+1],linewidth=4)
			if orientation == 'v':
				plt.plot(pts[i][j][:,0]-arr_n, pts[i][j][:,1], color=colorComp[lastComp+i+1],linewidth=4)
				plt.plot(xel[i][j]-arr_n, yelH[i][j], color=colorComp[lastComp+i+1],linewidth=4) 
				plt.plot(xelH[i][j]-arr_n, yel[i][j], color=colorComp[lastComp+i+1],linewidth=4)
	for j in xrange(0,len(xeltemp)): #newcomp
			arr_n = arrtemp[j]
			idx = int(indexComptemp[j])
			if orientation == 'h':
				plt.plot(ptstemp[j][:,0], ptstemp[j][:,1]-arr_n, color='red',linewidth=4)
				plt.plot(xeltemp[j], yelHtemp[j]-arr_n, color='red',linewidth=4) 
				plt.plot(xelHtemp[j], yelHtemp[j]-arr_n, color='red',linewidth=4)
			if orientation == 'v':
				plt.plot(ptstemp[j][:,0]-arr_n, ptstemp[j][:,1], color='red',linewidth=4)
				plt.plot(xeltemp[j]-arr_n, yelHtemp[j], color='red',linewidth=4) 
				plt.plot(xelHtemp[j]-arr_n, yelHtemp[j], color='red',linewidth=4)
		
def replace_comp(index,map_num,limplot_x1,limplot_x2,limplot_y1,limplot_y2,pts_arr,x_el_arr,x_elH_arr,
y_elH_arr,y_el_arr,realDAT,ext,arr,x,first_contour):
	print 'You selected a component in the ', map_num, 'map'
	print 'The component believed to be the same was selected automatically for the next maps'

	var_track = False
	while var_track == False:
		prompt = 'Are you happy with the automatical selection? (Y/n)'	
		var = str(raw_input(prompt))

		if 'Y' in var or 'y' in var or 'n' in var or len(var)==0:
			var_track = True	
	
	if var == 'n':
		changes = 'True'
		for i in xrange(0,len(index)):
			if (i+1) != map_num:
				print 'Do you want to change the component in map', i+1, '?'
				var_track = False
				while var_track == False:
					prompt = '(Y/n)'	
					var = str(raw_input(prompt))
	
					if 'Y' in var or 'y' in var or 'n' in var or len(var)==0:
						var_track = True	
				if 'Y' in var or 'y' in var or len(var)==0:
					print 'Select the new component'
					plt.ion()
					plt.figure(2)	
					plt.axis('scaled')
					plt.xlim(limplot_x1,limplot_x2)
					plt.ylim(limplot_y2,limplot_y1)	
    					plot_components(pts_arr,x_el_arr,x_elH_arr,y_elH_arr,y_el_arr,arr)
					plot_related_components(pts_arr,x_el_arr,x_elH_arr,y_elH_arr,y_el_arr,arr,index,'red')
					plot_maps(realDAT,ext,arr,first_contour)
					plt.xlabel('Right Ascension [mas]')

					#select the component
					param = ginput(2,0)	 
					x_c = float(param[1][0])
					near_comp = int(find_nearest(x[i],x_c)[1])
					index[i] = near_comp

					plt.close('all')
	else:
		changes = 'False'

	return index,changes


