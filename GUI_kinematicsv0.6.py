import numpy as np
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from scipy.optimize import curve_fit
import pylab as plt
import random
from pylab import *
from scipy import io
from scipy import stats
from scipy.optimize import leastsq
import scipy.optimize as optimization
import matplotlib.ticker as ticker
from math import *
import cmath as math
import pickle
import iminuit, probfit
import astropy.io.fits as pf
from searchNED import searchNEDnoGUI
from functions2 import read_map, saver
from functions2 import Annotate, convolve_difmap, convolve_difmap_modelfit, search_modelfit
from functions_conv import order_by_nu, conv_params, read_conv_params, check_map_params
from kinematics_functions import plot_components, plot_maps, plot_NiA_components, axis_limits#, searchNed#
from kinematics_functions import search_rms,take_header, read_modfile, x_y, order_by_date, find_nearest, shift_modelfit_difmap
from kinematics_functions import ellipse_axis_lines, plot_related_components, replace_comp,r_psi
from plot_components import get_ellipse_coords, ellipse_axis
import os,glob,re
import subprocess as sub
import urllib2

#if mpl.__version__ == '2.0.0':
#	mpl.style.use('classic')

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def plotting():
		color_comp = ['#F2A5EA','#A6E157','Cyan','MediumVioletRed','Gold','#A52A2A','#A25CF0','#240B3B','#393B0B','#FE2E64']
		#plotting not core components
		for k in xrange(0,len(KinematicsWindow.WINcompIndex)):
			components_index = KinematicsWindow.WINcompIndex[k]
			for i in xrange(0,len(components_index)):
				for j in xrange(0,len(KinematicsWindow.WINx_el_arr)):
					arr_n = needed_param.arr[i]
					idx = int(components_index[i])
					if KinematicsWindow.orientation == 'h':
						plt.plot(KinematicsWindow.WINpts_arr[i][idx][:,0], KinematicsWindow.WINpts_arr[i][idx][:,1]-arr_n, color=color_comp[k],linewidth=3)
      						plt.plot(KinematicsWindow.WINx_el_arr[i][idx], KinematicsWindow.WINy_elH_arr[i][idx]-arr_n, color=color_comp[k],linewidth=3) 	
       						plt.plot(KinematicsWindow.WINx_elH_arr[i][idx], KinematicsWindow.WINy_el_arr[i][idx]-arr_n, color=color_comp[k],linewidth=3)
					if KinematicsWindow.orientation == 'v':
						plt.plot(KinematicsWindow.WINpts_arr[i][idx][:,0]-arr_n, KinematicsWindow.WINpts_arr[i][idx][:,1], color=color_comp[k],linewidth=3)
      						plt.plot(KinematicsWindow.WINx_el_arr[i][idx]-arr_n, KinematicsWindow.WINy_elH_arr[i][idx], color=color_comp[k],linewidth=3) 	
       						plt.plot(KinematicsWindow.WINx_elH_arr[i][idx]-arr_n, KinematicsWindow.WINy_el_arr[i][idx], color=color_comp[k],linewidth=3)
			last_comp = k #for colors

		for i in xrange(0,len(KinematicsWindow.WINx_el_arrNiA)):
			arrn= KinematicsWindow.WINarrNiA[i]
			idxList = KinematicsWindow.WINcompAutomaticIndexNotInAll[i]
			for j in xrange(0,len(KinematicsWindow.WINx_el_arrNiA[i])):
				arr_n = arrn[j]
				idx = int(idxList[j])
				if KinematicsWindow.orientation == 'h':
					plt.plot(KinematicsWindow.WINpts_arrNiA[i][j][:,0], KinematicsWindow.WINpts_arrNiA[i][j][:,1]-arr_n, color=color_comp[last_comp+i+1],linewidth=4)
					plt.plot(KinematicsWindow.WINx_el_arrNiA[i][j], KinematicsWindow.WINy_elH_arrNiA[i][j]-arr_n, color=color_comp[last_comp+i+1],linewidth=4) 
					plt.plot(KinematicsWindow.WINx_elH_arrNiA[i][j], KinematicsWindow.WINy_el_arrNiA[i][j]-arr_n, color=color_comp[last_comp+i+1],linewidth=4)
				if KinematicsWindow.orientation == 'v':
					plt.plot(KinematicsWindow.WINpts_arrNiA[i][j][:,0]-arr_n, KinematicsWindow.WINpts_arrNiA[i][j][:,1], color=color_comp[last_comp+i+1],linewidth=4)
					plt.plot(KinematicsWindow.WINx_el_arrNiA[i][j]-arr_n, KinematicsWindow.WINy_elH_arrNiA[i][j], color=color_comp[last_comp+i+1],linewidth=4) 
					plt.plot(KinematicsWindow.WINx_elH_arrNiA[i][j]-arr_n, KinematicsWindow.WINy_el_arrNiA[i][j], color=color_comp[last_comp+i+1],linewidth=4)

		#plot new component
		if KinematicsWindow.NotInAll == 0:
			for i in xrange(0,len(KinematicsWindow.WINx_el_arr)):
				arr_n = needed_param.arr[i]
				idx = int(KinematicsWindow.WINcompAutomaticIndex[i])
				if KinematicsWindow.orientation == 'h':
					plt.plot(KinematicsWindow.WINpts_arr[i][idx][:,0], KinematicsWindow.WINpts_arr[i][idx][:,1]-arr_n, color='red',linewidth=4)
       					plt.plot(KinematicsWindow.WINx_el_arr[i][idx], KinematicsWindow.WINy_elH_arr[i][idx]-arr_n, color='red',linewidth=4) 
       					plt.plot(KinematicsWindow.WINx_elH_arr[i][idx], KinematicsWindow.WINy_el_arr[i][idx]-arr_n, color='red',linewidth=4)
				if KinematicsWindow.orientation == 'v':
					plt.plot(KinematicsWindow.WINpts_arr[i][idx][:,0]-arr_n, KinematicsWindow.WINpts_arr[i][idx][:,1], color='red',linewidth=4)
       					plt.plot(KinematicsWindow.WINx_el_arr[i][idx]-arr_n, KinematicsWindow.WINy_elH_arr[i][idx], color='red',linewidth=4) 
       					plt.plot(KinematicsWindow.WINx_elH_arr[i][idx]-arr_n, KinematicsWindow.WINy_el_arr[i][idx], color='red',linewidth=4)

		if KinematicsWindow.NotInAll == 1: #WORKING HERE
			for j in xrange(0,len(KinematicsWindow.WINx_el_arrNiAtemp)): #newcomp
					arr_n = KinematicsWindow.WINarrNiAtemp[j]
					idx = int(KinematicsWindow.WINcompAutomaticIndexNiAnewComp[j])
					if KinematicsWindow.orientation == 'h':
						plt.plot(KinematicsWindow.WINpts_arrNiAtemp[j][:,0], KinematicsWindow.WINpts_arrNiAtemp[j][:,1]-arr_n, color='red',linewidth=4)
       						plt.plot(KinematicsWindow.WINx_el_arrNiAtemp[j], KinematicsWindow.WINy_elH_arrNiAtemp[j]-arr_n, color='red',linewidth=4) 
       						plt.plot(KinematicsWindow.WINx_elH_arrNiAtemp[j], KinematicsWindow.WINy_el_arrNiAtemp[j]-arr_n, color='red',linewidth=4)
					if KinematicsWindow.orientation == 'v':
						plt.plot(KinematicsWindow.WINpts_arrNiAtemp[j][:,0]-arr_n, KinematicsWindow.WINpts_arrNiAtemp[j][:,1], color='red',linewidth=4)
       						plt.plot(KinematicsWindow.WINx_el_arrNiAtemp[j]-arr_n, KinematicsWindow.WINy_elH_arrNiAtemp[j], color='red',linewidth=4) 
       						plt.plot(KinematicsWindow.WINx_elH_arrNiAtemp[j]-arr_n, KinematicsWindow.WINy_el_arrNiAtemp[j], color='red',linewidth=4)


def line(t, m, c): 
	return m * t - c*m

def parabolic(t, m, c, acc,tmid): # define it to be parabolic or whatever you like
	return m * t - c*m + acc*(t - tmid)**2

class SearchNEDclass(QWidget):
	def __init__(self):
     	 	QWidget.__init__(self,None, Qt.WindowStaysOnTopHint)
	   	#self.layout = QGridLayout()
	   	self.layout = QGridLayout()
		#needed_param.source_name = '0836+710'
		self.sourceName = needed_param.source_name


		self.searchNED()
	

	def searchNED(self):
		if self.sourceName != '2013+370':
			if self.sourceName.find('+') > 0:
				splitedName = self.sourceName.split('+')
				name1 = splitedName[0]
				name2 = splitedName[1]
				response = urllib2.urlopen('https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname='+str(name1)+'%2B'+str(name2)+'&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES')
			if self.sourceName.find('-') > 0:
				splitedName = self.sourceName.split('-')
				name1 = splitedName[0]
				name2 = splitedName[1]
				response = urllib2.urlopen('https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname='+str(name1)+'-'+str(name2)+'&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES')

			html = response.read()

			DLpos = html.find('Luminosity Distance')
			if DLpos == -1:


				self.labelempty = QLabel()
				self.labelTEXT = QLabel()
				self.labelSOURCE = QLabel()
				self.labelTEXT2 = QLabel()
				self.labelTEXT.setText("Sorry the source was not found in NED by ")
				self.labelTEXT.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
				self.labelSOURCE.setText(self.sourceName)
				self.labelSOURCE.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
				self.labelTEXT2.setText("Please give a more common name or set DL and z manually")
				self.labelTEXT2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

				self.layout.addWidget(self.labelempty,1,0)
				self.layout.addWidget(self.labelempty,2,0)
				self.layout.addWidget(self.labelTEXT,1,1,1,3)
				self.layout.addWidget(self.labelSOURCE,1,4)
				self.layout.addWidget(self.labelTEXT2,2,1,1,5)

				self.labelSourceName = QLabel()
				self.labelSourceName.setText("Source name    : ")
				self.labelSourceName.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
				self.SourceName = QLineEdit()
				#self.SigmaCut.setValidator(QDoubleValidator())
				#self.SigmaCut.textChanged.connect(self.check_state)
				#self.SigmaCut.textChanged.emit(self.SigmaCut.text())
				self.SourceName.setFixedSize(100,25)

				self.labelDL = QLabel()
				self.labelDL.setText("DL                          : ")
				self.labelDL.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
				self.DL = QLineEdit()
				self.DL.setValidator(QDoubleValidator())
				self.DL.textChanged.connect(self.check_state)
				self.DL.textChanged.emit(self.DL.text())
				self.DL.setFixedSize(100,25)

				self.labelZ = QLabel()
				self.labelZ.setText("z                             : ")
				self.labelZ.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
				self.Z = QLineEdit()
				self.Z.setValidator(QDoubleValidator())
				self.Z.textChanged.connect(self.check_state)
				self.Z.textChanged.emit(self.Z.text())
				self.Z.setFixedSize(100,25)

				self.layout.addWidget(self.labelSourceName,3,1)
				self.layout.addWidget(self.SourceName,3,2)
				self.layout.addWidget(self.labelDL,4,1)
				self.layout.addWidget(self.DL,4,2)
				self.layout.addWidget(self.labelZ,5,1)
				self.layout.addWidget(self.Z,5,2)

				self.labelSourceName.setBuddy(self.SourceName)
				self.labelDL.setBuddy(self.DL)
				self.labelZ.setBuddy(self.Z)

				self.getSource = QPushButton("&Get")
				self.getSource.setFixedSize(100,25)
				self.getSource.clicked.connect(lambda: self.getParameters())

				self.layout.addWidget(self.getSource,6,2)

				self.setLayout(self.layout)
	
				DL = 0.
				z = 0.
				raDeg,decDeg = 0., 0.

			else:
				DLpos1st = DLpos+29
				afterDLpos = html[DLpos:]
				DLposEnd = DLpos+afterDLpos.find('Mpc')-2

				#print html[DLpos1st], html[DLposEnd]

				DL= ''
				for i in xrange (DLpos1st,DLposEnd+1):
					DL= DL + html[i]

				DL = float(DL)
				print(DL)

				z1 = html.find('HREF="#ObjNo1"')+13
				z2 = html[z1:]
				z3 = z1 + z2.find('&nbsp')
				z4 = html[z1:z3].split(' ')
				z5 = []
				for item in z4:
					if item:
						z5.append(item)

				z = float(z5[len(z5)-1])

				print 'DL =', DL, 'z=',z

				ra = z5[len(z5)-5]
				dec = z5[len(z5)-4]

				raSplit1 = ra.split('h')
				raSplit2 = raSplit1[1].split('m')
				raH = float(raSplit1[0])
				raM = float(raSplit2[0])
				raS = float(raSplit2[1].split('s')[0])	
				raHours = raH + raM/60 + raS/(60*60)
				raDeg = raHours*15

				decSplit1 = dec.split('d')
				decSplit2 = decSplit1[1].split('m')
				decD = float(decSplit1[0])
				decM = float(decSplit2[0])
				decS = float(decSplit2[1].split('s')[0])	
				if self.sourceName.find('+') > 0:
					decDeg = decD + decM/(60) + decS/(60*60)
				if self.sourceName.find('-') > 0:
					decDeg = decD - decM/(60) - decS/(60*60)

				needed_param.DL = DL
				needed_param.z = z

				#c = SkyCoord(ra=raHours*u.degree,dec=decDays*u.degree)
	
				print 'RA =', ra, 'DEC =', dec
				print 'RA =', raDeg, 'degrees, DEC =', decDeg, 'degrees' 

				self.closeIt()
		else:
			DL = 5489 #mpc #for 2013+370
			z = 0.859
			raDeg = 0.
			decDeg = 0.

			print 'DL =', DL, 'z=',z

			needed_param.DL = DL
			needed_param.z = z

	def check_state(self,*args,**kwargs):
		sender = self.sender()
		validator = sender.validator()
		state = validator.validate(sender.text(),0)[0]
		if state == QValidator.Acceptable:
			color = '#c4df9b' #green
		elif state == QValidator.Intermediate:
			color = '#fff79a' #yellow
		else:
			color = '#f6989d' #red
		sender.setStyleSheet('QLineEdit { background-color: %s }' %color)

	def getParameters(self):
		if self.SourceName.text() != '':
			self.sourceName = str(self.SourceName.text())

			self.searchNED()

    			self.closeIt()
		else:
			needed_param.DL = float(self.DL.text())
			needed_param.z = float(self.Z.text())

			print 'DL =', needed_param.DL, 'z=',needed_param.z
      			self.closeIt()

	def closeIt(self): 
		self.close()

class SelectEpochs(QWidget):
	def __init__(self):
     	 	QWidget.__init__(self,None, Qt.WindowStaysOnTopHint)
	   	layout = QGridLayout()
		#layout.addStretch(1)
		#layout.addLayout(hbox)
     	  	
		self.labelempty = QLabel()
		self.labelTEXT = QLabel()
		self.labelTEXT2 = QLabel()
		self.checks = []
          	
		self.labelTEXT.setText("Select the epochs in which you want to modify the selected component")
		self.labelTEXT.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

		layout.addWidget(self.labelempty,3,0)
		layout.addWidget(self.labelTEXT,0,1,1,5)

        	for i in xrange(0,len(needed_param.date)):
            		c = QCheckBox('%s' %(needed_param.date[i],))
	    		c.setFixedSize(100,25)
            		layout.addWidget(c,i+1,1)
            		self.checks.append(c)
	   
        	self.SelectButton = QPushButton("&Select \n manually")
		self.SelectButton.setAutoDefault(True)

		for i in xrange(0,len(needed_param.date)+5):
			layout.addWidget(self.labelempty, i, 0)	
		if KinematicsWindow.NotInAll == 1:
		        self.SelectButtonAuto = QPushButton("&Select \n automatically")
			self.SelectButtonAuto.setAutoDefault(True)
			layout.addWidget(self.SelectButtonAuto, len(needed_param.date)+4,1,2,1)

		layout.addWidget(self.SelectButton, len(needed_param.date)+2,1,2,1)

		plt.ioff()
		self.plot = plt.figure()		
		self.canvas = FigureCanvas(self.plot)
		self.toolbar = NavigationToolbar(self.canvas,self)
		self.Plot_All()
		layout.addWidget(self.canvas,1,2,len(needed_param.date)+5,len(needed_param.date)+5)
		layout.addWidget(self.toolbar,len(needed_param.date)+7,2,1,len(needed_param.date)+5)

		self.setLayout(layout)

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

	def Plot_All(self):
		color_comp = ['#F2A5EA','#A6E157','Cyan','MediumVioletRed','Gold','#A52A2A','#A25CF0','#240B3B','#393B0B','#FE2E64']
    		for i in xrange(0,len(KinematicsWindow.WINx_el_arr)):
			arr_n = needed_param.arr[i]
			ax = self.plot.add_subplot(111)
			for j in xrange(0,len(KinematicsWindow.WINx_el_arr[i])):
				if KinematicsWindow.orientation == 'h':	
        				ax.plot(KinematicsWindow.WINpts_arr[i][j][:,0], KinematicsWindow.WINpts_arr[i][j][:,1]-arr_n, color='blue',linewidth=3)
        				ax.plot(KinematicsWindow.WINx_el_arr[i][j], KinematicsWindow.WINy_elH_arr[i][j]-arr_n, color='blue',linewidth=3) 
        				ax.plot(KinematicsWindow.WINx_elH_arr[i][j], KinematicsWindow.WINy_el_arr[i][j]-arr_n, color='blue',linewidth=3)
				if KinematicsWindow.orientation == 'v':	
        				ax.plot(KinematicsWindow.WINpts_arr[i][j][:,0]-arr_n, KinematicsWindow.WINpts_arr[i][j][:,1], color='blue',linewidth=3)
        				ax.plot(KinematicsWindow.WINx_el_arr[i][j]-arr_n, KinematicsWindow.WINy_elH_arr[i][j], color='blue',linewidth=3) 
        				ax.plot(KinematicsWindow.WINx_elH_arr[i][j]-arr_n, KinematicsWindow.WINy_el_arr[i][j], color='blue',linewidth=3)
		for i in xrange(0,len(needed_param.realDAT)):
			arr_n = needed_param.arr[i]
			ax = self.plot.add_subplot(111)
			levels = needed_param.first_contour[i]*needed_param.realDAT[i].std()*np.array([-1., 1., 1.41,2.,2.83,4.,5.66,8.,11.3,16.,
		                                22.6,32.,45.3,64.,90.5,128.,181.,256.,362.,512.,
		                                724.,1020.,1450.,2050.])
			if KinematicsWindow.orientation == 'h':
				ext = [needed_param.ext[0],needed_param.ext[1],needed_param.ext[2]-arr_n,needed_param.ext[3]-arr_n]
			if KinematicsWindow.orientation == 'v':
				ext = [needed_param.ext[0]-arr_n,needed_param.ext[1]-arr_n,needed_param.ext[2],needed_param.ext[3]]
			cset = ax.contour(needed_param.realDAT[i], levels, inline=1,
			                  colors=['grey']
			                  ,extent=ext, aspect=1.0
			                  )
		plotting()
		ax.axis('scaled')
		ax.set_xlim(KinematicsWindow.WINlimits[0],KinematicsWindow.WINlimits[1])
		ax.set_ylim(KinematicsWindow.WINlimits[3],KinematicsWindow.WINlimits[2])
		if KinematicsWindow.orientation == 'h': 
			plt.xlabel('Right Ascension [mas]')
			ax = plt.gca()
			ax.set_yticks(needed_param.arr)
			ax.set_yticklabels(needed_param.date[::-1])
		if KinematicsWindow.orientation == 'v': 
			plt.ylabel('Declination [mas]')
			ax = plt.gca()
			ax.set_xticks(needed_param.arr)
			ax.set_xticklabels(needed_param.date[::-1])
		self.canvas.draw()

class popupAutomaticSelection(QWidget):
	def __init__(self):
     	 	QWidget.__init__(self)
	   	layout = QGridLayout()
		#layout.addStretch(1)
		#layout.addLayout(hbox)
     	  	
		self.labelempty = QLabel()
		self.labelTEXT = QLabel()
		self.labelTEXT2 = QLabel()
          	
		self.labelTEXT.setText("The components for the rest of the maps have been selected automatically")
		self.labelTEXT.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		self.labelTEXT2.setText("Are you happy with the automatic selection?")
		self.labelTEXT2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

		layout.addWidget(self.labelempty,3,0)
		layout.addWidget(self.labelTEXT,0,1,1,5)
		layout.addWidget(self.labelTEXT2,1,1,1,5)
	   
        	self.YesButton = QPushButton("&Yes")
        	self.NoButton = QPushButton("&No")
		self.YesButton.setAutoDefault(True)
		self.NoButton.setAutoDefault(True)
		
		layout.addWidget(self.YesButton, 3,2)
		layout.addWidget(self.NoButton, 3,3)

		#for the plot in the window 
		plt.ioff()
		self.plot = plt.figure()		
		self.canvas = FigureCanvas(self.plot)
		self.toolbar = NavigationToolbar(self.canvas,self)
		self.Plot_All()
		layout.addWidget(self.canvas,4,0,9,9)
		layout.addWidget(self.toolbar,15,0,1,9)

		#self.canvas = FigureCanvas(KinematicsWindow.plot)
       		#layout.addWidget(self.canvas)

		self.setLayout(layout)

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

	def Plot_All(self):
		color_comp = ['#F2A5EA','#A6E157','Cyan','MediumVioletRed','Gold','#A52A2A','#A25CF0','#240B3B','#393B0B','#FE2E64']
    		for i in xrange(0,len(KinematicsWindow.WINx_el_arr)):
			arr_n = needed_param.arr[i]
			ax = self.plot.add_subplot(111)
			for j in xrange(0,len(KinematicsWindow.WINx_el_arr[i])):
				if KinematicsWindow.orientation == 'h':	
        				ax.plot(KinematicsWindow.WINpts_arr[i][j][:,0], KinematicsWindow.WINpts_arr[i][j][:,1]-arr_n, color='blue',linewidth=3)
        				ax.plot(KinematicsWindow.WINx_el_arr[i][j], KinematicsWindow.WINy_elH_arr[i][j]-arr_n, color='blue',linewidth=3) 
        				ax.plot(KinematicsWindow.WINx_elH_arr[i][j], KinematicsWindow.WINy_el_arr[i][j]-arr_n, color='blue',linewidth=3)
				if KinematicsWindow.orientation == 'v':	
        				ax.plot(KinematicsWindow.WINpts_arr[i][j][:,0]-arr_n, KinematicsWindow.WINpts_arr[i][j][:,1], color='blue',linewidth=3)
        				ax.plot(KinematicsWindow.WINx_el_arr[i][j]-arr_n, KinematicsWindow.WINy_elH_arr[i][j], color='blue',linewidth=3) 
        				ax.plot(KinematicsWindow.WINx_elH_arr[i][j]-arr_n, KinematicsWindow.WINy_el_arr[i][j], color='blue',linewidth=3)
		for i in xrange(0,len(needed_param.realDAT)):
			arr_n = needed_param.arr[i]
			ax = self.plot.add_subplot(111)
			levels = needed_param.first_contour[i]*needed_param.realDAT[i].std()*np.array([-1., 1., 1.41,2.,2.83,4.,5.66,8.,11.3,16.,
		                                22.6,32.,45.3,64.,90.5,128.,181.,256.,362.,512.,
		                                724.,1020.,1450.,2050.])
			if KinematicsWindow.orientation == 'h':
				ext = [needed_param.ext[0],needed_param.ext[1],needed_param.ext[2]-arr_n,needed_param.ext[3]-arr_n]
			if KinematicsWindow.orientation == 'v':
				ext = [needed_param.ext[0]-arr_n,needed_param.ext[1]-arr_n,needed_param.ext[2],needed_param.ext[3]]
			cset = ax.contour(needed_param.realDAT[i], levels, inline=1,
			                  colors=['grey']
			                  ,extent=ext, aspect=1.0
			                  )
		plotting()
		ax.axis('scaled')
		ax.set_xlim(KinematicsWindow.WINlimits[0],KinematicsWindow.WINlimits[1])
		ax.set_ylim(KinematicsWindow.WINlimits[3],KinematicsWindow.WINlimits[2])
		if KinematicsWindow.orientation == 'h': 
			plt.xlabel('Right Ascension [mas]')
			ax = plt.gca()
			ax.set_yticks(needed_param.arr)
			ax.set_yticklabels(needed_param.date[::-1])
		if KinematicsWindow.orientation == 'v': 
			plt.ylabel('Declination [mas]')
			ax = plt.gca()
			ax.set_xticks(needed_param.arr)
			ax.set_xticklabels(needed_param.date[::-1])
		self.canvas.draw()


class popupLegendColors(QWidget):
	def __init__(self):
     	 	super(popupLegendColors, self).__init__()
	   	self.layout = QGridLayout()
		#layout.addStretch(1)
		#layout.addLayout(hbox)

		self.colorsRect = []
		self.colorsPU = []
		self.labelsPU = []
		self.zOrderPU = []
     	  	
		self.CompName = QLabel()
		self.ColorLabel = QLabel()
		self.CompName.setText("Name")
		self.CompName.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		self.CompName.setFixedSize(100,25)
		self.ColorLabel.setText("Color")
		self.ColorLabel.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		self.ColorLabel.setFixedSize(100,25)

		self.layout.addWidget(self.CompName,0,1)
		self.layout.addWidget(self.ColorLabel,0,2)

		for i in xrange(0,len(KinematicsWindow.Labels)):
		    c = QLineEdit()
		    c.setText(KinematicsWindow.Labels[i])
		    c.setFixedSize(100,25)
		    self.labelsPU.append(c)
		    c2 = QLineEdit()
		    c2.setText(KinematicsWindow.colorsFinal[i])
		    c2.setFixedSize(100,25)
		    self.colorsPU.append(c2)
		    self.layout.addWidget(c,i+1,1)
		    self.layout.addWidget(c2,i+1,2)
		    c3 = QLabel()
		    c3.setFixedSize(50,20)
		    c3.setStyleSheet('QLabel {background-color:'+str(KinematicsWindow.colorsFinal[i])+'; border-color: 20px solid black}')
		    self.layout.addWidget(c3,i+1,3)
		    c4 = QLineEdit()
		    c4.setText('%d' % KinematicsWindow.zOrderPlot[i])
		    c4.setFixedSize(50,20)
		    self.zOrderPU.append(c4)
		    self.layout.addWidget(c4,i+1,4)

        	self.PickButton = QPushButton("&Set Labels \n and colors")
		self.PickButton.setAutoDefault(True)
		
		self.layout.addWidget(self.PickButton, len(KinematicsWindow.Labels)+1,2,2,1)

		#for the plot in the window 
		plt.ioff()
		self.plot = plt.figure()		
		self.canvas = FigureCanvas(self.plot)
		self.toolbar = NavigationToolbar(self.canvas,self)
		self.Plot_All()
		self.layout.addWidget(self.canvas,1,5,9,6)
		self.layout.addWidget(self.toolbar,10,5,1,6)

		self.setLayout(self.layout)

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

	def Plot_All(self):

		"""		if np.all(np.asarray(self.xComp[i])<0) or np.all(np.asarray(self.yComp[i])<0) :
			yaxis = np.asarray(self.rComp[i]).copy()
		else:
			yaxis =- np.asarray(self.rComp[i]).copy()"""

		ax = self.plot.add_subplot(111)
		for i in xrange(len(KinematicsWindow.WINr[0])):
			if np.all(KinematicsWindow.WINx[0][i]<0) or np.all(KinematicsWindow.WINy[0][i]<0) :
				yaxis = KinematicsWindow.WINr[0][i]
			else:
				yaxis = -KinematicsWindow.WINr[0][i]		

			plt.errorbar(needed_param.epoch,yaxis,yerr=KinematicsWindow.WINrerr[0][i],color=KinematicsWindow.colorsFinal[i],fmt='o')
			#plt.errorbar(needed_param.epoch,KinematicsWindow.WINr[0][i],color=KinematicsWindow.colorsFinal[i],marker='o',linestyle='None')
			last = i 
		for i in xrange(len(KinematicsWindow.WINrNiA[0])):
			if np.all(np.asarray(KinematicsWindow.WINxNiA[0][i])<0) or np.all(np.asarray(KinematicsWindow.WINyNiA[0][i])<0) :
				yaxis = np.asarray(KinematicsWindow.WINrNiA[0][i])
			else:
				yaxis = -np.asarray(KinematicsWindow.WINrNiA[0][i])	
			plt.errorbar(KinematicsWindow.WINcompAutomaticIndexNotInAllepoch[i],yaxis,yerr=KinematicsWindow.WINrNiAerr[0][i],color=KinematicsWindow.colorsFinal[last+i+1],fmt='o')
		ax2 = plt.gca()
		ax2.get_xaxis().get_major_formatter().set_useOffset(False)
		ax2.minorticks_on()
		ax2.tick_params('both', length=10, width=2, which='major')
		ax2.tick_params('both', length=5, width=1, which='minor')
		#ax.set_ylim(bottom=ymin)
		plt.xlabel('Time [years]')
		plt.ylabel('r [mas]')
		self.canvas.draw()


	def check_state(self,*args,**kwargs):
		sender = self.sender()
		validator = sender.validator()
		state = validator.validate(sender.text(),0)[0]
		if state == QValidator.Acceptable:
			color = '#c4df9b' #green
		elif state == QValidator.Intermediate:
			color = '#fff79a' #yellow
		else:
			color = '#f6989d' #red
		sender.setStyleSheet('QLineEdit { background-color: %s }' %color)


class KinematicsWindow(QWidget): 

    NotInAll = 0
    recheck = 0
    automatic = 0

    WINx = np.asarray([])
    WINy = np.asarray([])

    #parameters for plotting the ellipse read from the needed_param class
    WINpts_arr = np.asarray([])
    WINx_el_arr,WINy_el_arr = np.asarray([]),np.asarray([])
    WINx_elH_arr,WINy_elH_arr = np.asarray([]),np.asarray([])

    WINcompIndex = np.asarray([])

    WINcompAutomaticIndex = []

    #for components that are not in all epochs
    WINpts_arrNiA = []
    WINx_el_arrNiA,WINy_el_arrNiA = [],[]
    WINx_elH_arrNiA,WINy_elH_arrNiA = [],[]
    WINarrNiA = []

    WINcompAutomaticIndexNotInAll = []

    WINcompAutomaticIndexNotInAllepoch = []


    #for components that are not in all epochs temporal	
    WINcompAutomaticIndexNiAnewComp = []
    WINcompAutomaticIndexNiAepochnewComp = []

    WINpts_arrNiAtemp = np.asarray([])
    WINx_el_arrNiAtemp,WINy_el_arrNiAtemp = np.asarray([]),np.asarray([])
    WINx_elH_arrNiAtemp,WINy_elH_arrNiAtemp = np.asarray([]),np.asarray([])
    WINarrNiAtemp = np.asarray([])
    WINr = np.asarray([])
    WINrNiA = np.asarray([])
    WINrerr = np.asarray([])
    WINrNiAerr = np.asarray([])
    WINx = np.asarray([])
    WINxNiA = np.asarray([])
    WINy = np.asarray([])
    WINyNiA = np.asarray([])

    WINlimits = []
    orientation = ' '

    Labels = []
    colorsFinal = []
    zOrderPlot = []


    def __init__(self,*args):
        QWidget.__init__(self)

        self.layout = QGridLayout()
	#layout.addStretch(1)
	#layout.addLayout(hbox)

	self.Labels = []
        self.ColorsFinal = []

	#initializing parameters

	print needed_param.DL, needed_param.z
        self.date1 = []
	self.date2 = []

	self.fits1 = 'l'
	self.fits2 = 'r'
	self.freq1 = 0.
	self.freq2 = 0.
	self.files_chosen = []
	self.models_chosen = [] 
	self.shifted_files = []

	#beam,mapsize and cellsize of all the maps, read from the needed_param class
	self.bmaj_files = needed_param.bmaj
	self.bmin_files = needed_param.bmin
	self.bpa_files = needed_param.bpa
	self.beam_files = needed_param.beam
	self.cells_files = needed_param.pixelsize
	self.mapsize_file = needed_param.mapsize

	#components parameters read from the needed_param class
	self.r, self.errr =   np.asarray(needed_param.r).copy(), np.asarray(needed_param.errr).copy()
	self.psi, self.errpsi =  np.asarray(needed_param.psi).copy(), np.asarray(needed_param.errpsi) .copy()
	self.size, self.errsize = np.asarray(needed_param.size).copy(), np.asarray(needed_param.errsize).copy()
	self.flux, self.errflux = np.asarray(needed_param.flux).copy(), np.asarray(needed_param.errflux).copy()
	self.tb, self.errtb = np.asarray(needed_param.tb).copy(), np.asarray(needed_param.errtb).copy()

	self.x, self.errx = np.asarray(needed_param.x).copy(), np.asarray(needed_param.errx).copy()
	self.y, self.erry = np.asarray(needed_param.y).copy(), np.asarray(needed_param.erry).copy()

	#parameters for plotting the ellipse read from the needed_param class
	self.pts_arr = np.asarray(needed_param.pts_arr).copy()
	self.x_el_arr,self.y_el_arr = np.asarray(needed_param.x_el_arr).copy(),np.asarray(needed_param.y_el_arr).copy()
	self.x_elH_arr,self.y_elH_arr = np.asarray(needed_param.x_elH_arr).copy(),np.asarray(needed_param.y_elH_arr).copy()

	#asigning to class parameters parameters that i need in the popup window
    	KinematicsWindow.WINx = self.x
    	KinematicsWindow.WINy = self.y
    	KinematicsWindow.WINpts_arr = self.pts_arr
    	KinematicsWindow.WINx_el_arr,KinematicsWindow.WINy_el_arr = self.x_el_arr,self.y_el_arr
    	KinematicsWindow.WINx_elH_arr,KinematicsWindow.WINy_elH_arr = self.x_elH_arr,self.y_elH_arr

	#orientation of the source
	self.orientation = 'h'
	KinematicsWindow.orientation = self.orientation

	#component parameters after shifting
	#to initialize them, they take the value of the initial arrays
	#this will change later in the program 
	#only r,psi, x and y change (the errors remain the same as well)
	self.rShifted = []
	self.errrShifted = []
	self.psiShifted = []
	self.errpsiShifted = []

	self.xNew = np.asarray(needed_param.x).copy()
	self.yNew = np.asarray(needed_param.y).copy()

	#arrays containing the shift values
	self.shiftX, self.shiftY = np.array([0.]*len(needed_param.x)), np.array([0.]*len(needed_param.x)) 

	self.rms1 = 0.
	self.rms2 = 0.
	self.ext = []

	#edges of all the maps, read from the needed_param class
	self.limplot_x1 = needed_param.ext[0]
	self.limplot_x2 = needed_param.ext[1]
	self.limplot_y1 = needed_param.ext[2]
	self.limplot_y2 = needed_param.ext[3]

	#array containing the position of the core for each epoch
	self.index_c = np.array([0.]*len(needed_param.x))

	#array containing the position of the component to be identified for each epoch
	self.indexIdentifying = np.array([0.]*len(needed_param.x))

	#because the files will be shifted later with respect of the core position
	#the names fo the files that will contain the shifted images are created
	self.shift_files = []
	for i in xrange(0,len(needed_param.fits)):
		self.shift_files.append('SHIFT/'+str('%1.2f' % (needed_param.epoch[i]))+'shifted.fits')

	#colors for the components
	self.colors = ['#F2A5EA','#A6E157','Cyan','MediumVioletRed','Gold','#A52A2A','#A25CF0','#240B3B','#393B0B','#FE2E64']

	#total number of subplots
	self.ax = []
	for i in xrange(0,len(needed_param.date)):
		self.ax.append('ax'+str(i+1))

	#lists with the ordered components
	self.rComp, self.errrComp = [], []
	self.psiComp, self.psierrComp = [], []
	self.sizeComp, self.sizeerrComp = [], []
	self.fluxComp, self.fluxerrComp = [], []
	self.tbComp, self.tberrComp = [], []
	self.xComp, self.xerrComp = [], []
	self.yComp, self.yerrComp = [], []

	self.rCompNiA, self.errrCompNiA = [], []
	self.psiCompNiA, self.psierrCompNiA = [], []
	self.sizeCompNiA, self.sizeerrCompNiA = [], []
	self.fluxCompNiA, self.fluxerrCompNiA = [], []
	self.tbCompNiA, self.tberrCompNiA = [], []
	self.xCompNiA, self.xerrCompNiA = [], []
	self.yCompNiA, self.yerrCompNiA = [], []

	#arrays for NiA components before saying that I agree with the selection
	self.rCompNiAtemp, self.errrCompNiAtemp = [], []
	self.psiCompNiAtemp, self.psierrCompNiAtemp = [], []
	self.sizeCompNiAtemp, self.sizeerrCompNiAtemp = [], []
	self.fluxCompNiAtemp, self.fluxerrCompNiAtemp = [], []
	self.tbCompNiAtemp, self.tberrCompNiAtemp = [], []
	self.xCompNiAtemp, self.xerrCompNiAtemp = [], []
	self.yCompNiAtemp, self.yerrCompNiAtemp= [], []

	#list with the fits results
	self.LineValuesX = []
	self.ParabolicValuesX = []
	self.LineValuesY = []
	self.ParabolicValuesY = []

	self.errLineValuesX = []
	self.errParabolicValuesX = []
	self.errLineValuesY = []
	self.errParabolicValuesY = []

	#list with the fits results
	self.LineValuesXnia = []
	self.ParabolicValuesXnia = []
	self.LineValuesYnia = []
	self.ParabolicValuesYnia = []

	self.errLineValuesXnia = []
	self.errParabolicValuesXnia = []
	self.errLineValuesYnia = []
	self.errParabolicValuesYnia = []

	#store the fit values for each component with best chi2
	self.FitValuesX = []
	self.FitValuesY = []
	self.errFitValuesX = []
	self.errFitValuesY = []

	#fits with the velocity results
	self.muLine, self.betaLine = [], []
	self.muerrLine, self.betaerrLine = [], []
	self.tejjLine, self.tejjLineerr = [], []

	self.mu, self.beta = [], []
	self.muerr, self.betaerr = [], []

	self.load = False
	self.changeColors = False

        for i in xrange(0,len(needed_param.date)):
		c1 = QLabel() 
		c1.setText('date %d:' %(i+1,))
		c1.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		c1.setFixedSize(50,25)
		c2 = QLabel() 
		c2.setText('%s' %(needed_param.date[i],))
		c2.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
		self.layout.addWidget(c1,i+1,1)
		self.layout.addWidget(c2,i+1,2,1,3)
		self.date1.append(c1)
		self.date2.append(c2)

	self.labelempty = QLabel()
	self.labelempty.setFixedSize(25,25)
	self.labelempty2 = QLabel()


	for i in xrange(0,10):
		self.layout.addWidget(self.labelempty, i, 0)	
		self.layout.addWidget(self.labelempty, i, 4)
		self.layout.addWidget(self.labelempty, i, 5)	
		self.layout.addWidget(self.labelempty, i, 6)	


   	self.labelSeparation = QLabel()
   	self.labelSeparation.setText("Separation:")
	self.Separation = QLineEdit()
	self.Separation.setValidator(QDoubleValidator())
	self.Separation.textChanged.connect(self.check_state)
	self.Separation.textChanged.emit(self.Separation.text())
	self.Separation.setFixedSize(100,25)
	self.Separation.setText('%1.1f' % (20.0))
	self.layout.addWidget(self.labelSeparation, 1, 7)
	self.layout.addWidget(self.Separation, 1, 8)

   	self.labelSigmaCut = QLabel()
   	self.labelSigmaCut.setText("Sigma Cut:")
	self.SigmaCut = QLineEdit()
	self.SigmaCut.setValidator(QDoubleValidator())
	self.SigmaCut.textChanged.connect(self.check_state)
	self.SigmaCut.textChanged.emit(self.SigmaCut.text())
	self.SigmaCut.setFixedSize(100,25)
	self.SigmaCut.setText('%1.1f' % (1.0))
	self.layout.addWidget(self.labelSigmaCut, 2, 7)
	self.layout.addWidget(self.SigmaCut, 2, 8)


   	self.labelDL = QLabel()
   	self.labelDL2 = QLabel()
   	self.labelDL.setText("DL =")
	self.labelDL.setFixedSize(50,25)
   	self.labelDL2.setText(str(needed_param.DL)+" Mpc")
	self.labelDL.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelDL2.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelDL.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
	self.labelDL2.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
	self.layout.addWidget(self.labelDL, 1, 10)
	self.layout.addWidget(self.labelDL2, 1, 11)

   	self.labelz = QLabel()
   	self.labelz2 = QLabel()
   	self.labelz.setText("z    =")
	self.labelz.setFixedSize(50,25)
   	self.labelz2.setText(str(needed_param.z))
	self.labelz.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelz2.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelz.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
	self.labelz2.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
	self.layout.addWidget(self.labelz, 2, 10)
	self.layout.addWidget(self.labelz2, 2, 11)
	
   	self.labelsource = QLabel()
   	self.labelsource2 = QLabel()
   	self.labelsource.setText("source: ")
   	self.labelsource2.setText(str(needed_param.source_name))
	self.labelsource.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelsource2.setStyleSheet('QLabel {color: magenta } QLabel {font: Bold }')
	self.layout.addWidget(self.labelsource, 0, 0)
	self.layout.addWidget(self.labelsource2, 0, 1)

	#self.progressBar = QProgressBar(self)
        #self.progressBar.setRange(0,100)

        self.SelectRegionButton = QPushButton("&Plot maps")
        self.SelectRegionButton.clicked.connect(lambda: self.SelectRegion())
	self.SelectRegionButton.setFixedSize(340,25)
	self.SelectRegionButton.setAutoDefault(True)

        self.SelectCoreButton = QPushButton("&Select Core")
        self.SelectCoreButton.clicked.connect(lambda: self.SelectCore())
	self.SelectCoreButton.setFixedSize(340,25)
	self.SelectCoreButton.setAutoDefault(True)

        self.SelectCompInAllButton = QPushButton("&Select Component in all epochs")
        self.SelectCompInAllButton.clicked.connect(lambda: self.SelectCompInAll())
	self.SelectCompInAllButton.setFixedSize(340,25)
	self.SelectCompInAllButton.setAutoDefault(True)

        self.SelectCompInSomeButton = QPushButton("&Select Component with missing epochs")
        self.SelectCompInSomeButton.clicked.connect(lambda: self.SelectCompInSome())
	self.SelectCompInSomeButton.setFixedSize(340,25)
	self.SelectCompInSomeButton.setAutoDefault(True)

        self.ReplaceCompButton = QPushButton("&Replace Component")
        self.ReplaceCompButton.clicked.connect(lambda: self.ReplaceComp())
	self.ReplaceCompButton.setFixedSize(340,25)
	self.ReplaceCompButton.setAutoDefault(True)

        self.SetAllButton = QPushButton("&Set All Components")
        self.SetAllButton.clicked.connect(lambda: self.SetAll(False))
	self.SetAllButton.setFixedSize(340,25)
	self.SetAllButton.setAutoDefault(True)

        self.FittingButton = QPushButton("&Fitting")
        self.FittingButton.clicked.connect(lambda: self.Fitting())
	self.FittingButton.setFixedSize(340,25)
	self.FittingButton.setAutoDefault(True)

        self.LoadButton = QPushButton("&Load existing identification")
        self.LoadButton.clicked.connect(lambda: self.Load())
	self.LoadButton.setFixedSize(200,25)
	self.LoadButton.setAutoDefault(True)

        self.ChangeColorsButton = QPushButton("&Change Label and Colors")
        self.ChangeColorsButton.clicked.connect(lambda: self.SetAll(True))
	self.ChangeColorsButton.setFixedSize(200,25)
	self.ChangeColorsButton.setAutoDefault(True)

	
	self.layout.addWidget(self.SelectRegionButton, 4, 7,1,5)
	self.layout.addWidget(self.SelectCoreButton, 5, 7,1,5)
	self.layout.addWidget(self.SelectCompInAllButton, 6, 7,1,5)
	self.layout.addWidget(self.SelectCompInSomeButton, 7, 7,1,5)
	self.layout.addWidget(self.ReplaceCompButton, 8, 7,1,5)
	self.layout.addWidget(self.SetAllButton, 9, 7,1,5)
	self.layout.addWidget(self.FittingButton, 10, 7,1,5)
	self.layout.addWidget(self.LoadButton,4, 13,1,3)
	self.layout.addWidget(self.ChangeColorsButton,5, 13,1,3)

	#layout.addWidget(self.progressBar,17,1,1,len(freq))"""

        self.setLayout(self.layout)

	#put the window in the center of the desktop
	qr = self.frameGeometry()
	cp = QDesktopWidget().availableGeometry().center()
	qr.moveCenter(cp)
	self.move(qr.topLeft())

	self.setWindowTitle("Kinematics")

    #function to check if the arguments given in the text boxes are fine
    #if the value is not the kind of parameter wanted, 
    #for example, if the box requires a double and you give an integrer, it get red
    #if the input value is still the kind of parameter wanted, but outside a range, it gets yellow
    #if the input value is fine, it gets green
    def check_state(self,*args,**kwargs):
		sender = self.sender()
		validator = sender.validator()
		state = validator.validate(sender.text(),0)[0]
		if state == QValidator.Acceptable:
			color = '#c4df9b' #green
		elif state == QValidator.Intermediate:
			color = '#fff79a' #yellow
		else:
			color = '#f6989d' #red
		sender.setStyleSheet('QLineEdit { background-color: %s }' %color)


    def SelectRegion(self):

	sigma_cut = float(self.SigmaCut.text())
	separation = float(self.Separation.text())

	needed_param.first_contour = []	
	for i in xrange(0,len(needed_param.rms_maps)):
		needed_param.first_contour.append(sigma_cut*needed_param.rms_maps[i])

	if len(needed_param.x)%2==0:
		lims = needed_param.lenght/2*separation
		needed_param.arr = np.linspace(-lims,lims,needed_param.lenght)
	else:
		lims = trunc(needed_param.lenght/2)*separation
		needed_param.arr = np.linspace(-lims,lims,needed_param.lenght)

	plt.figure(1)
	plt.axis('scaled')
	plt.xlim(200,-200) 
	plt.ylim(needed_param.arr[0]-160,needed_param.arr[len(needed_param.arr)-1]+160)
    	plot_components(self.pts_arr,self.x_el_arr,self.x_elH_arr,self.y_elH_arr,self.y_el_arr,needed_param.arr,self.orientation)
	plot_maps(needed_param.realDAT,needed_param.ext,needed_param.arr,needed_param.first_contour,self.orientation)
	if self.orientation == 'h': 
		plt.xlabel('Right Ascension [mas]')
		ax = plt.gca()
		ax.set_yticks(needed_param.arr)
		ax.set_yticklabels(needed_param.date[::-1])
	if self.orientation == 'v': 
		plt.ylabel('Declination [mas]')
		ax = plt.gca()
		ax.set_xticks(needed_param.arr)
		ax.set_xticklabels(needed_param.date[::-1])

	a = Annotate()
	plt.show()

	#self.layout.addWidget(self.canvas, 3,8, 8,8)
	#self.layout.addWidget(self.toolbar,11,8,8,8)
		
	[self.limplot_x1,self.limplot_x2,self.limplot_y1,self.limplot_y2] = a()
	KinematicsWindow.WINlimits = [self.limplot_x1,self.limplot_x2,self.limplot_y1,self.limplot_y2]


    def SelectCore(self):

	#[self.limplot_x1,self.limplot_x2,self.limplot_y1,self.limplot_y2] = [10.7780612245, -48.9158163265, 55.1020408163, -50.5102040816] #[30.4464285714, -11.0841836735, 70.7142857143, -79.693877551][10.7780612245, -48.9158163265, 55.1020408163, -50.5102040816]
	#KinematicsWindow.WINlimits = [self.limplot_x1,self.limplot_x2,self.limplot_y1,self.limplot_y2]
	plt.ion()
	plt.figure(2)	
	plt.axis('scaled')
	plt.xlim(self.limplot_x1,self.limplot_x2)
	plt.ylim(self.limplot_y2,self.limplot_y1)	
    	plot_components(self.pts_arr,self.x_el_arr,self.x_elH_arr,self.y_elH_arr,self.y_el_arr,needed_param.arr,self.orientation)
	plot_maps(needed_param.realDAT,needed_param.ext,needed_param.arr,needed_param.first_contour,self.orientation)
	if self.orientation == 'h': 
		plt.xlabel('Right Ascension [mas]')
		ax = plt.gca()
		ax.set_yticks(needed_param.arr)
		ax.set_yticklabels(needed_param.date[::-1])
	if self.orientation == 'v': 
		plt.ylabel('Declination [mas]')
		ax = plt.gca()
		ax.set_xticks(needed_param.arr)
		ax.set_xticklabels(needed_param.date[::-1])

	#select the core component
	param = ginput(2,0) 
	reversed_arr = needed_param.arr[::-1]
	if self.orientation == 'h':
		x_c = float(param[1][0])
		ind = int(find_nearest(reversed_arr,float(param[1][1]))[1])
		map_num = ind+1
	if self.orientation == 'v':
		y_c = float(param[1][1])
		ind = int(find_nearest(reversed_arr,float(param[1][0]))[1])
		map_num = ind+1

	plt.close('all')

	#find the nearest component to the one selected in all maps
	
	old_comp = -100.
	for i in xrange(0,len(needed_param.x)):
		if self.orientation == 'h':
			near_comp = int(find_nearest(self.x[i],x_c)[1])
		if self.orientation == 'v':
			near_comp = int(find_nearest(self.y[i],y_c)[1])
		#while  old_comp < near_comp:
		#	near_comp = int(find_nearest(x[i],x_c+1.)[1])
		self.index_c[i] = int(near_comp)

	#putting the core component list in the list that has all the components
	KinematicsWindow.WINcompIndex = np.asarray([self.index_c])

	#putting the core component list in the list that has the component just asigned automatically
	#this will change each time with each new component asigned
	KinematicsWindow.WINcompAutomaticIndex = self.index_c

	self.PopUpSelection()

    #only for the core selection
    def PopUpSelection(self):
	self.wi = popupAutomaticSelection()
	self.wi.show()
	self.wi.YesButton.clicked.connect(lambda: self.GoodCoreSelection())
        self.wi.NoButton.clicked.connect(lambda: self.BadCoreSelection())

    def GoodCoreSelection(self):
	"""	self.rShifted = []
	self.errrShifted = []
	self.psiShifted = []
	self.errpsiShifted = []

	self.xNew = np.asarray(needed_param.x).copy()
	self.yNew = np.asarray(needed_param.y).copy()

	#arrays containing the shift values
	self.shiftX, self.shiftY = np.array([0.]*len(x)), np.array([0.]*len(x)) """

	for i in xrange(0,len(self.x)):
		shiftx = self.x[i][int(self.index_c[i])]
		shifty = self.y[i][int(self.index_c[i])]
		self.shiftX[i] = shiftx
		self.shiftY[i] = shifty
		print 'The shift is: RA = ', shiftx, 'mas, DEC = ', shifty, 'mas' 
		#convolve_difmap([files[i]],[models[i]],bmaj,bmin,bpa,shiftx,shifty,mapsize,pixelsize,[shift_files[i]])
		convolve_difmap_modelfit([needed_param.files[i]],[needed_param.models[i]],[needed_param.modelfit[i]],self.bmaj_files,self.bmin_files,self.bpa_files,-shiftx,-shifty,self.mapsize_file,self.cells_files,[self.shift_files[i]])
		search_modelfit()
		#remove log file
		os.system('rm difmap.log*\n')
		#getting the shifted values for r and psi for all the components and epochs
		mod_parameters = read_modfile(['newModelFit.txt'],self.beam_files,needed_param.bmajReal,needed_param.bminReal, needed_param.z,needed_param.freq,needed_param.rms_maps/1000.,'Done')
		self.rShifted.append(mod_parameters[0][0]) 
		self.errrShifted = self.errr.copy() 
		self.psiShifted.append(mod_parameters[2][0])
		self.errpsiShifted = self.errpsi.copy()
   
		#remove mod file
		os.system('rm newModelFit.txt*\n')

		#shifting x and y values for all the components in all epochs
		x_l = []
		y_l = []
		for j in xrange(0,len(self.x[i])):
			x_l.append(self.x[i][j]-shiftx)	
			y_l.append(self.y[i][j]-shifty)    

		self.xNew[i] = np.asarray(x_l)
		self.yNew[i] = np.asarray(y_l)

	#recalculating the axis and outer circle for the plot of the ellipses components
	ellipse_plot = ellipse_axis_lines(self.xNew,self.yNew,self.size)
	self.pts_arr,self.pt_arr = np.asarray(ellipse_plot[0]), np.asarray(ellipse_plot[1])
	self.x_el_arr,self.y_el_arr = np.asarray(ellipse_plot[2]), np.asarray(ellipse_plot[3])
	self.x_elH_arr,self.y_elH_arr = np.asarray(ellipse_plot[4]), np.asarray(ellipse_plot[5]) 

	#giving these parameters also to the window class parameters
	KinematicsWindow.WINx = self.xNew
	KinematicsWindow.WINy = self.yNew
	KinematicsWindow.WINpts_arr = self.pts_arr
	KinematicsWindow.WINx_el_arr,KinematicsWindow.WINy_el_arr = self.x_el_arr,self.y_el_arr
	KinematicsWindow.WINx_elH_arr,KinematicsWindow.WINy_elH_arr = self.x_elH_arr,self.y_elH_arr

	KinematicsWindow.WINcompIndex = np.asarray([self.index_c])

	#KinematicsWindow.Labels.append('core')

	self.wi.close()
	
    def BadCoreSelection(self):
	isCore = 1
	self.wi2 = SelectEpochs()
	self.wi2.show()
	self.wi2.SelectButton.clicked.connect(lambda: self.ReselectComps(self.wi2.checks,self.x,self.y,isCore,0))

    def ReselectComps(self,checkBoxes,x,y,isCore,automatic):

	if KinematicsWindow.NotInAll == 1:
		self.wi2.close()

	plt.ion()
	plt.figure(4)	
	plt.axis('scaled')
	plt.xlim(self.limplot_x1,self.limplot_x2)
	plt.ylim(self.limplot_y2,self.limplot_y1)
    	plot_components(self.pts_arr,self.x_el_arr,self.x_elH_arr,self.y_elH_arr,self.y_el_arr,needed_param.arr,self.orientation)
	plotting()
	plot_maps(needed_param.realDAT,needed_param.ext,needed_param.arr,needed_param.first_contour,self.orientation)
	if self.orientation == 'h': 
		plt.xlabel('Right Ascension [mas]')
		ax = plt.gca()
		ax.set_yticks(needed_param.arr)
		ax.set_yticklabels(needed_param.date[::-1])
	if self.orientation == 'v': 
		plt.ylabel('Declination [mas]')
		ax = plt.gca()
		ax.set_xticks(needed_param.arr)
		ax.set_xticklabels(needed_param.date[::-1])

	epochSelected = []
	for i in xrange(len(checkBoxes)):
		if checkBoxes[i].isChecked():
			epochSelected.append(needed_param.epoch[i])

	if automatic == 0:
		param = ginput(2*len(epochSelected),0) 
		map_numC = []

		j = -1

	if automatic == 1:
		param = ginput(2,0) 
		reversed_arr = needed_param.arr[::-1]
		if self.orientation == 'h':
			x_comp = float(param[1][0])
			ind = int(find_nearest(reversed_arr,float(param[1][1]))[1])
		if self.orientation == 'v':
			y_comp = float(param[1][1])
			ind = int(find_nearest(reversed_arr,float(param[1][0]))[1])
		map_num = ind+1

		index_comp_l = []
		old_comp = -100.

		components_index = KinematicsWindow.WINcompIndex

	k = -1
	j = -1
	lists = []

	n = []
	for i in xrange(len(checkBoxes)):
		if checkBoxes[i].isChecked():
			n.append(i)

	tempCompIndex = np.asarray([0.]*len(n))
	tempCompEpoch = np.asarray([0.]*len(n))
	temparr = np.asarray([0.]*len(n))
	pts1 = np.empty(self.pts_arr[0][0].shape)#,np.zeros(len(self.pts_arr[0][0][1]))])
	xel = np.empty(self.x_el_arr[0][0].shape)
	tempPts = np.empty(len(n), dtype=object)
	tempXel = np.empty(len(n), dtype=object)
	tempYel = np.empty(len(n), dtype=object)
	tempXelH = np.empty(len(n), dtype=object)
	tempYelH = np.empty(len(n), dtype=object)
	for i in xrange(0,len(n)):
		tempPts[i] = pts1
		tempXel[i] = xel
		tempYel[i] = xel
		tempXelH[i] = xel
		tempYelH[i] = xel

	if KinematicsWindow.recheck == 0:
		rNiA = np.asarray([0.]*len(n))
		errrNiA = np.asarray([0.]*len(n))
		psiNiA = np.asarray([0.]*len(n))
		errpsiNiA = np.asarray([0.]*len(n))
		sizeNiA = np.asarray([0.]*len(n))
		errsizeNiA = np.asarray([0.]*len(n))
		fluxNiA = np.asarray([0.]*len(n))
		errfluxNiA = np.asarray([0.]*len(n))
		tbNiA = np.asarray([0.]*len(n))
		errtbNiA = np.asarray([0.]*len(n))
		xNiA = np.asarray([0.]*len(n))
		errxNiA = np.asarray([0.]*len(n))
		yNiA = np.asarray([0.]*len(n))
		erryNiA = np.asarray([0.]*len(n))

	if KinematicsWindow.recheck == 1:
		rNiA = np.asarray(self.rCompNiAtemp)
		errrNiA = np.asarray(self.errrCompNiAtemp)
		psiNiA = np.asarray(self.psiCompNiAtemp)
		errpsiNiA = np.asarray(self.psierrCompNiAtemp)
		sizeNiA = np.asarray(self.sizeCompNiAtemp)
		errsizeNiA = np.asarray(self.sizeerrCompNiAtemp)
		fluxNiA = np.asarray(self.fluxCompNiAtemp)
		errfluxNiA = np.asarray(self.fluxerrCompNiAtemp)
		tbNiA = np.asarray(self.tbCompNiAtemp)
		errtbNiA = np.asarray(self.tberrCompNiAtemp)
		xNiA = np.asarray(self.xCompNiAtemp)
		errxNiA = np.asarray(self.xerrCompNiAtemp)
		yNiA = np.asarray(self.yCompNiAtemp)
		erryNiA = np.asarray(self.yerrCompNiAtemp)

	for i in xrange(len(checkBoxes)):
		if checkBoxes[i].isChecked():
			if automatic == 0:
				j= j + 2
				k = k +1	
				if self.orientation == 'h':
					x_Comp = float(param[j][0])
					ind = int(find_nearest(x[i],x_Comp)[1])
				if self.orientation == 'v':
					y_Comp = float(param[j][1])
					ind = int(find_nearest(y[i],y_Comp)[1])
				if isCore == 1:
					self.index_c[i] = ind

			if automatic == 1:
				k = k + 1
				if self.orientation == 'h':
					near_comp = int(find_nearest(x[i],x_comp)[1])
					x_pos_comp = x[i][int(find_nearest(x[i],x_comp)[1])]
				if self.orientation == 'v':
					near_comp = int(find_nearest(y[i],y_comp)[1])
					y_pos_comp = y[i][int(find_nearest(y[i],y_comp)[1])]
				steps = np.linspace(0,100,1000,endpoint=True)
				temp = []
				for j in xrange(len(components_index)):
					temp.append(KinematicsWindow.WINcompIndex[j][i]) #defined by [component][epoch]
				tracker = -1
				while near_comp in temp:
					tracker = tracker +1
					if self.orientation == 'h':
						if x_pos_comp < 0:
							near_comp = int(find_nearest(self.xNew[i],x_comp-steps[tracker])[1])
						else:
							near_comp = int(find_nearest(self.xNew[i],x_comp+steps[tracker])[1])
					if self.orientation == 'v':
						if y_pos_comp < 0:
							near_comp = int(find_nearest(self.yNew[i],y_comp-steps[tracker])[1])
						else:
							near_comp = int(find_nearest(self.yNew[i],y_comp+steps[tracker])[1])
				ind = near_comp

			if KinematicsWindow.NotInAll == 1:
				if KinematicsWindow.recheck == 0:
					tempCompIndex[k] = ind
					tempCompEpoch[k] = needed_param.epoch[i]
					temparr[k] = needed_param.arr[i]
					tempPts[k] = np.asarray(self.pts_arr[i][ind])
					tempXel[k] = np.asarray(self.x_el_arr[i][ind])
					tempYel[k] = np.asarray(self.y_el_arr[i][ind])
					tempXelH[k] = np.asarray(self.x_elH_arr[i][ind])
					tempYelH[k] = np.asarray(self.y_elH_arr[i][ind])

					#put a small temp parameter and in the second elif put that list equal to the ones that are here right now

					rNiA[k] = self.rShifted[i][ind]
					errrNiA[k] = self.errrShifted[i][ind]
					psiNiA[k] = self.psiShifted[i][ind]
					errpsiNiA[k] = self.errpsiShifted[i][ind]
					sizeNiA[k] = self.size[i][ind]
					errsizeNiA[k] = self.errsize[i][ind]
					fluxNiA[k] = self.flux[i][ind]
					errfluxNiA[k] = self.errflux[i][ind]
					tbNiA[k] = self.tb[i][ind]
					errtbNiA[k] = self.errtb[i][ind]
					xNiA[k] = self.xNew[i][ind]
					errxNiA[k] = self.errx[i][ind]
					yNiA[k] = self.yNew[i][ind]
					erryNiA[k] = self.erry[i][ind]

				if KinematicsWindow.recheck == 1:

					replace_index = np.where(KinematicsWindow.WINcompAutomaticIndexNiAepochnewComp==needed_param.epoch[i])[0]

					KinematicsWindow.WINcompAutomaticIndexNiAnewComp[replace_index[0]] = ind
					KinematicsWindow.WINcompAutomaticIndexNiAepochnewComp[replace_index[0]] = needed_param.epoch[i]
					KinematicsWindow.WINarrNiAtemp[replace_index[0]] = needed_param.arr[i]
			    		KinematicsWindow.WINpts_arrNiAtemp[replace_index[0]] = np.asarray(self.pts_arr[i][ind])
			    		KinematicsWindow.WINx_el_arrNiAtemp[replace_index[0]] = np.asarray(self.x_el_arr[i][ind])
					KinematicsWindow.WINy_el_arrNiAtemp[replace_index[0]] = np.asarray(self.y_el_arr[i][ind])
			    		KinematicsWindow.WINx_elH_arrNiAtemp[replace_index[0]] = np.asarray(self.x_elH_arr[i][ind])
					KinematicsWindow.WINy_elH_arrNiAtemp[replace_index[0]] = np.asarray(self.y_elH_arr[i][ind])

					rNiA[replace_index] = self.rShifted[i][ind]
					errrNiA[replace_index] = self.errrShifted[i][ind]
					psiNiA[replace_index] = self.psiShifted[i][ind]
					errpsiNiA[replace_index] = self.errpsiShifted[i][ind]
					sizeNiA[replace_index] = self.size[i][ind]
					errsizeNiA[replace_index] = self.errsize[i][ind]
					fluxNiA[replace_index] = self.flux[i][ind]
					errfluxNiA[replace_index] = self.errflux[i][ind]
					tbNiA[replace_index] = self.tb[i][ind]
					errtbNiA[replace_index] = self.errtb[i][ind]
					xNiA[replace_index] = self.xNew[i][ind]
					errxNiA[replace_index] = self.errx[i][ind]
					yNiA[replace_index] = self.yNew[i][ind] 
					erryNiA[replace_index] = self.erry[i][ind]


			else:
				self.indexIdentifying[i] = ind
		#else:
		#	if KinematicsWindow.NotInAll == 1:
		#		lists.append('nan')

	if isCore == 1:
		#putting the core component list in the list that has all the components
		KinematicsWindow.WINcompIndex = np.asarray([self.index_c])

		#putting the core component list in the list that has the component just asigned automatically
		#this will change each time with each new component asigned
		KinematicsWindow.WINcompAutomaticIndex = self.index_c
	elif KinematicsWindow.NotInAll == 1:
		if KinematicsWindow.recheck == 0:
			KinematicsWindow.WINcompAutomaticIndexNiAnewComp = tempCompIndex
			KinematicsWindow.WINcompAutomaticIndexNiAepochnewComp = tempCompEpoch
			KinematicsWindow.WINarrNiAtemp = temparr
	    		KinematicsWindow.WINpts_arrNiAtemp = tempPts
	    		KinematicsWindow.WINx_el_arrNiAtemp = tempXel
			KinematicsWindow.WINy_el_arrNiAtemp = tempYel
	    		KinematicsWindow.WINx_elH_arrNiAtemp = tempXelH
			KinematicsWindow.WINy_elH_arrNiAtemp = tempYelH

		self.rCompNiAtemp, self.errrCompNiAtemp = rNiA.tolist(), errrNiA.tolist() 
		self.psiCompNiAtemp, self.psierrCompNiAtemp = psiNiA.tolist(), errpsiNiA.tolist() 
		self.sizeCompNiAtemp, self.sizeerrCompNiAtemp = sizeNiA.tolist(), errsizeNiA.tolist() 
		self.fluxCompNiAtemp, self.fluxerrCompNiAtemp = fluxNiA.tolist(), errfluxNiA.tolist() 
		self.tbCompNiAtemp, self.tberrCompNiAtemp = tbNiA.tolist(), errtbNiA.tolist() 
		self.xCompNiAtemp, self.xerrCompNiAtemp = xNiA.tolist(), errxNiA.tolist() 
		self.yCompNiAtemp, self.yerrCompNiAtemp = yNiA.tolist(), erryNiA.tolist() 

			
	else:
		KinematicsWindow.WINcompAutomaticIndex = self.indexIdentifying


	plt.close('all')
	self.wi.close()
	self.wi2.close()

	if isCore == 1:
		self.PopUpSelection()

	else:
		if KinematicsWindow.NotInAll == 1:
			self.PopUpSelectionOther()
			KinematicsWindow.recheck = 1
		else:		
			self.PopUpSelectionOther()

    def SelectCompInAll(self):
	KinematicsWindow.NotInAll = 0

	plt.ion()
	plt.figure(4)	
	plt.axis('scaled')
	plt.xlim(self.limplot_x1,self.limplot_x2)
	plt.ylim(self.limplot_y2,self.limplot_y1)
    	plot_components(self.pts_arr,self.x_el_arr,self.x_elH_arr,self.y_elH_arr,self.y_el_arr,needed_param.arr,self.orientation)
	components_index = KinematicsWindow.WINcompIndex	
	plotting()
	plot_maps(needed_param.realDAT,needed_param.ext,needed_param.arr,needed_param.first_contour,self.orientation)
	if self.orientation == 'h': 
		plt.xlabel('Right Ascension [mas]')
		ax = plt.gca()
		ax.set_yticks(needed_param.arr)
		ax.set_yticklabels(needed_param.date[::-1])
	if self.orientation == 'v': 
		plt.ylabel('Declination [mas]')
		ax = plt.gca()
		ax.set_xticks(needed_param.arr)
		ax.set_xticklabels(needed_param.date[::-1])

	param = ginput(2,0) 
	reversed_arr = needed_param.arr[::-1]
	if self.orientation == 'h':
		x_comp = float(param[1][0])
		ind = int(find_nearest(reversed_arr,float(param[1][1]))[1])
	if self.orientation == 'v':
		y_comp = float(param[1][1])
		ind = int(find_nearest(reversed_arr,float(param[1][0]))[1])
	map_num = ind+1

	index_comp_l = []
	old_comp = -100.
	for i in xrange(0,len(self.xNew)):
		if self.orientation == 'h':
			near_comp = int(find_nearest(self.xNew[i],x_comp)[1])
			x_pos_comp = self.xNew[i][int(find_nearest(self.xNew[i],x_comp)[1])]
		if self.orientation == 'v':
			near_comp = int(find_nearest(self.yNew[i],y_comp)[1])
			y_pos_comp = self.yNew[i][int(find_nearest(self.yNew[i],y_comp)[1])]
		steps = np.linspace(0,100,1000,endpoint=True)
		temp = []
		#for k in xrange(0,len(KinematicsWindow.WINcompIndex)):
		#	components_index = KinematicsWindow.WINcompIndex[k]
		for j in xrange(len(components_index)):
			temp.append(components_index[j][i]) 
		tracker = -1
		while near_comp in temp:
			tracker = tracker +1
			if self.orientation == 'h':
				if x_pos_comp < 0:
					near_comp = int(find_nearest(self.xNew[i],x_comp-steps[tracker])[1])
				else:
					near_comp = int(find_nearest(self.xNew[i],x_comp+steps[tracker])[1])
			if self.orientation == 'v':
				if y_pos_comp < 0:
					near_comp = int(find_nearest(self.yNew[i],y_comp-steps[tracker])[1])
				else:
					near_comp = int(find_nearest(self.yNew[i],y_comp+steps[tracker])[1])
		index_comp_l.append(near_comp)
		self.indexIdentifying = np.asarray(index_comp_l)

	KinematicsWindow.WINcompAutomaticIndex = self.indexIdentifying 

	plt.close('all')
	self.PopUpSelectionOther()

    def PopUpSelectionOther(self):
	self.wi = popupAutomaticSelection()
	self.wi.show()
	self.wi.YesButton.clicked.connect(lambda: self.GoodCompSelection())
        self.wi.NoButton.clicked.connect(lambda: self.BadCompSelection())

    def GoodCompSelection(self):
	if KinematicsWindow.NotInAll == 0:
		KinematicsWindow.WINcompIndex = np.vstack([KinematicsWindow.WINcompIndex,np.asarray([self.indexIdentifying])])
		#KinematicsWindow.Labels.append('c'+str(len(KinematicsWindow.Labels)))
	if KinematicsWindow.NotInAll == 1:
		KinematicsWindow.WINcompAutomaticIndexNotInAll.append(KinematicsWindow.WINcompAutomaticIndexNiAnewComp)
		KinematicsWindow.WINcompAutomaticIndexNotInAllepoch.append(KinematicsWindow.WINcompAutomaticIndexNiAepochnewComp)
		KinematicsWindow.WINarrNiA.append(KinematicsWindow.WINarrNiAtemp)
    		KinematicsWindow.WINpts_arrNiA.append(KinematicsWindow.WINpts_arrNiAtemp)
    		KinematicsWindow.WINx_el_arrNiA.append(KinematicsWindow.WINx_el_arrNiAtemp)
		KinematicsWindow.WINy_el_arrNiA.append(KinematicsWindow.WINy_el_arrNiAtemp)
    		KinematicsWindow.WINx_elH_arrNiA.append(KinematicsWindow.WINx_elH_arrNiAtemp)
		KinematicsWindow.WINy_elH_arrNiA.append(KinematicsWindow.WINy_elH_arrNiAtemp) 

		self.rCompNiA.append(self.rCompNiAtemp)
		self.errrCompNiA.append(self.errrCompNiAtemp)
		self.psiCompNiA.append(self.psiCompNiAtemp)
		self.psierrCompNiA.append(self.psierrCompNiAtemp)
		self.sizeCompNiA.append(self.sizeCompNiAtemp)
		self.sizeerrCompNiA.append(self.sizeerrCompNiAtemp)
		self.fluxCompNiA.append(self.fluxCompNiAtemp)
		self.fluxerrCompNiA.append(self.fluxerrCompNiAtemp)
		self.tbCompNiA.append(self.tbCompNiAtemp)
		self.tberrCompNiA.append(self.tberrCompNiAtemp)
		self.xCompNiA.append(self.xCompNiAtemp)
		self.xerrCompNiA.append(self.xerrCompNiAtemp)
		self.yCompNiA.append(self.yCompNiAtemp)
		self.yerrCompNiA.append(self.yerrCompNiAtemp)

		#KinematicsWindow.Labels.append('c'+str(len(KinematicsWindow.Labels)))

		KinematicsWindow.WINcompAutomaticIndexNiAnewComp = []
		KinematicsWindow.WINcompAutomaticIndexNiAepochnewComp = []
		KinematicsWindow.WINarrNiAtemp = []
    		KinematicsWindow.WINpts_arrNiAtemp = []
    		KinematicsWindow.WINx_el_arrNiAtemp = []
		KinematicsWindow.WINy_el_arrNiAtemp = []
    		KinematicsWindow.WINx_elH_arrNiAtemp = []
		KinematicsWindow.WINy_elH_arrNiAtemp = []

	self.wi.close()

    def BadCompSelection(self):
	isCore = 0
	self.wi2 = SelectEpochs()
	self.wi2.show()
	self.wi2.SelectButton.clicked.connect(lambda: self.ReselectComps(self.wi2.checks,self.xNew,self.yNew,isCore,0))

    def SelectCompInSome(self):
	KinematicsWindow.NotInAll = 1
	KinematicsWindow.recheck = 0
	isCore = 0
	self.wi2 = SelectEpochs()
	self.wi2.show()
	self.wi2.SelectButton.clicked.connect(lambda: self.ReselectComps(self.wi2.checks,self.xNew,self.yNew,isCore,0))
	self.wi2.SelectButtonAuto.clicked.connect(lambda: self.AutomaticSelectionInSome(self.wi2.checks))
	#self.wi = popupAutomaticSelection()
	#self.wi.show()

    def AutomaticSelectionInSome(self,checksBox): 

	self.wi2.close()

	KinematicsWindow.NotInAll = 1
	#onlyselectsone, the first one

	
	isCore = 0
	#KinematicsWindow.automatic = 1

	self.ReselectComps(self.wi2.checks,self.xNew,self.yNew,isCore,1)


    def ReplaceComp(self):
	print 'hi5' 

    def PopUpLegend(self):
	self.wi = popupLegendColors()
	self.wi.show()
	self.wi.PickButton.clicked.connect(lambda: self.Picking())


    def Picking(self):


	KinematicsWindow.Labels, KinematicsWindow.colorsFinal, KinematicsWindow.zOrderPlot = [], [], []

	for i in xrange(0,len(self.wi.colorsPU)):
		KinematicsWindow.Labels.append(str(self.wi.labelsPU[i].text()))
		KinematicsWindow.colorsFinal.append(str(self.wi.colorsPU[i].text()))
		KinematicsWindow.zOrderPlot.append(int(self.wi.zOrderPU[i].text()))

	#print KinematicsWindow.Labels
	#print KinematicsWindow.colorsFinal

	self.wi.close()

    def SetAll(self,changeColor):
	"""rComp, errrComp = [], []
	psiComp, psierrComp = [], []
	sizeComp, sizeerrComp = [], []
	xComp, xerrComp = [], []
	yComp, yerrComp = [], []"""


	if self.load == False:
		#defining temporal list to append to the big ordered list
		r_c, errr_c = [], []
		psi_c, errpsi_c = [], []
		size_c, errsize_c = [], []
		flux_c, errflux_c = [], []
		tb_c, errtb_c = [], []

		x_comp, y_comp = [], []
		xerr_comp, yerr_comp = [],[]

		self.rComp, self.errrComp = [], []
		self.psiComp, self.psierrComp = [], []
		self.sizeComp, self.sizeerrComp = [], []
		self.fluxComp, self.fluxerrComp = [], []
		self.tbComp, self.tberrComp = [], []
		self.xComp, self.xerrComp = [], []
		self.yComp, self.yerrComp = [], []

		components_index = KinematicsWindow.WINcompIndex
		#components_index2 = 
		for i in xrange(len(components_index)):
			for j in xrange(len(components_index[i])):
				r_c.append(self.rShifted[j][int(components_index[i][j])])
				errr_c.append(self.errrShifted[j][int(components_index[i][j])])
				psi_c.append(self.psiShifted[j][int(components_index[i][j])])
				errpsi_c.append(self.errpsiShifted[j][int(components_index[i][j])])
				size_c.append(self.size[j][int(components_index[i][j])])
				errsize_c.append(self.errsize[j][int(components_index[i][j])])
				flux_c.append(self.flux[j][int(components_index[i][j])])
				errflux_c.append(self.errflux[j][int(components_index[i][j])])
				tb_c.append(self.tb[j][int(components_index[i][j])])
				errtb_c.append(self.errtb[j][int(components_index[i][j])])
				x_comp.append(self.xNew[j][int(components_index[i][j])])
				y_comp.append(self.yNew[j][int(components_index[i][j])])
				xerr_comp.append(self.errx[j][int(components_index[i][j])])
				yerr_comp.append(self.erry[j][int(components_index[i][j])])
			

	#	self.rCompNiA, self.errrCompNiA = [], []
	#	self.psiCompNiA, self.psierrCompNiA = [], []
	#	self.sizeCompNiA, self.sizeerrCompNiA = [], []
	#	self.xCompNiA, self.xerrCompNiA = [], []
	#	self.yCompNiA, self.yerrCompNiA = [], []

			self.rComp.append(r_c)
			self.errrComp.append(errr_c)
			self.psiComp.append(psi_c)
			self.psierrComp.append(errpsi_c)
			self.sizeComp.append(size_c)
			self.sizeerrComp.append(errsize_c)
			self.fluxComp.append(flux_c)
			self.fluxerrComp.append(errflux_c)
			self.tbComp.append(tb_c)
			self.tberrComp.append(errtb_c)
			self.xComp.append(x_comp)
			self.xerrComp.append(xerr_comp)
			self.yComp.append(y_comp)
			self.yerrComp.append(yerr_comp)

			r_c, errr_c = [], []
			psi_c, errpsi_c = [], []
			size_c, errsize_c = [], []
			flux_c, errflux_c = [], []
			tb_c, errtb_c = [], []
			x_comp, xerr_comp = [], []
			y_comp, yerr_comp = [], []


		"""for i in xrange(0,len(self.rCompNiA)):
			self.rComp.append(self.rCompNiA[i])
			self.errrComp.append(self.errrCompNiA[i])
			self.psiComp.append(self.psiCompNiA[i])
			self.psierrComp.append(self.psierrCompNiA[i])
			self.sizeComp.append(self.sizeCompNiA[i])
			self.sizeerrComp.append(self.sizeerrCompNiA[i])
			self.xComp.append(self.xCompNiA[i])
			self.xerrComp.append(self.xerrCompNiA[i])
			self.yComp.append(self.yCompNiA[i])
			self.yerrComp.append(self.yerrCompNiA[i])"""

		plt.close('all')

	if self.load == True:
		components_index = self.rComp

	KinematicsWindow.WINr = np.asarray([self.rComp]).copy()
	KinematicsWindow.WINrNiA = np.asarray([self.rCompNiA]).copy()
	KinematicsWindow.WINrerr = np.asarray([self.errrComp]).copy()
	KinematicsWindow.WINrNiAerr = np.asarray([self.errrCompNiA]).copy()
	KinematicsWindow.WINx = np.asarray([self.xComp]).copy()
	KinematicsWindow.WINxNiA = np.asarray([self.xCompNiA]).copy()
	KinematicsWindow.WINy = np.asarray([self.yComp]).copy()
	KinematicsWindow.WINyNiA = np.asarray([self.yCompNiA]).copy()

	n = len(self.rComp)+len(self.rCompNiA)

	#print n

	if self.changeColors == False:
		KinematicsWindow.Labels = []

		KinematicsWindow.Labels.append('core')

		for i in xrange(1,n):
			KinematicsWindow.Labels.append('c'+str(len(KinematicsWindow.Labels)))

		KinematicsWindow.colorsFinal = self.colors[0:n]

		KinematicsWindow.zOrderPlot = np.arange(1,len(KinematicsWindow.colorsFinal)+1)

	if changeColor == True:
		self.PopUpLegend()
		self.changeColors = True

	#if self.load == False:
	filenames = []
	for i in xrange(0,len(self.fluxComp)):
		final_txt1 = [needed_param.epoch,self.fluxComp[i],self.fluxerrComp[i],self.rComp[i], self.errrComp[i],self.psiComp[i], self.psierrComp[i],self.sizeComp[i],self.sizeerrComp[i],self.xComp[i],self.xerrComp[i],self.yComp[i],self.yerrComp[i],self.tbComp[i],self.tberrComp[i]]
		header = np.array([['#Epoch (yr) Flux (Jy/beam) errFlux (Jy/beam) r (mas) errR (mas) psi (degrees) errPsi (degrees) size (mas) errSize (mas) RA (mas) errRA (mas) DEC (mas) errDEC (mas) Tb errTb \n #Component '+str(KinematicsWindow.Labels[i])]])
		final_txt=np.swapaxes(final_txt1, 0,1)
		file_name = 'ParametersFiles/Parameters'+str(i)+'.txt'
		filenames.append(file_name)	
		saver(file_name, header, final_txt, FORMAT='%1.4f')
		lasti = i

	for i in xrange(0,len(self.fluxCompNiA)):
		final_txt1 = [KinematicsWindow.WINcompAutomaticIndexNotInAllepoch[i],self.fluxCompNiA[i],self.fluxerrCompNiA[i],self.rCompNiA[i], self.errrCompNiA[i],self.psiCompNiA[i], self.psierrCompNiA[i],self.sizeCompNiA[i],self.sizeerrCompNiA[i],self.xCompNiA[i],self.xerrCompNiA[i],self.yCompNiA[i],self.yerrCompNiA[i],self.tbCompNiA[i],self.tberrCompNiA[i]]
		header = np.array([['#Epoch (yr) Flux (Jy/beam) errFlux (Jy/beam) r (mas) errR (mas) psi (degrees) errPsi (degrees) size (mas) errSize (mas) RA (mas) errRA (mas) DEC (mas) errDEC (mas) Tb errTb \n #Component '+str(KinematicsWindow.Labels[len(self.rComp)+i])]])
		final_txt=np.swapaxes(final_txt1, 0,1)
		file_name = 'ParametersFiles/NiAParameters'+str(lasti+i+1)+'.txt'
		filenames.append(file_name)	
		saver(file_name, header, final_txt, FORMAT='%1.4f')


	with open('All_parameters.txt', 'w') as outfile:
	    for fname in filenames:
		with open(fname) as infile:
		    for line in infile:
			outfile.write(line)

	self.Plott()

    def Plott(self):

	components_index = self.xComp
	plotMarginX = (np.max(needed_param.epoch)-np.min(needed_param.epoch))*1/3.8
	LabelsTempA = np.asarray(KinematicsWindow.Labels).copy()
	LabelsTemp = LabelsTempA.tolist()
	ColorsTempA = np.asarray(KinematicsWindow.colorsFinal).copy()
	ColorsTemp = LabelsTempA.tolist()
	ZTempA = np.asarray(KinematicsWindow.zOrderPlot).copy()
	ZTemp = LabelsTempA.tolist()
	zorderA = []

	#interpolation 
	#components_index = KinematicsWindow.WINcompIndex
	#for i in xrange(0,len(components_index)):
	#	x = linspace(self.xComp[i][0],self.xComp[i][len(self.xComp[i])-1],100,endpoint=False)
	#	y = self.yComp[i][0] + (x-self.xComp[i][0])*(self.yComp[i][len(self.yComp[i])-1]-self.yComp[i][0])/(self.xComp[i][len(self.xComp[i])-1]-self.xComp[i][0])
	#y = y1 +(x-x1)*(y2-y1)/(x2-x1)


	plt.figure(10)
	for i in xrange(len(components_index)):
		zorderA.append(plt.errorbar(needed_param.epoch,self.xComp[i],yerr=self.xerrComp[i],color=KinematicsWindow.colorsFinal[i],fmt='o'))#,label=KinematicsWindow.Labels[i])#,zorder=KinematicsWindow.zOrderPlot[i])
		last_comp = i
	for i in xrange(len(self.rCompNiA)):
		zorderA.append(plt.errorbar(KinematicsWindow.WINcompAutomaticIndexNotInAllepoch[i],self.xCompNiA[i],yerr=self.xerrCompNiA[i],color=KinematicsWindow.colorsFinal[last_comp+i+1],fmt='o'))#,label=KinematicsWindow.Labels[last_comp+i+1])#,zorder=KinematicsWindow.zOrderPlot[last_comp+i+1])
	#print ZTempA
	for i in range(len(ZTempA)):
		for j in range(1,len(ZTemp)-i):
			if ZTempA[j-1] > ZTempA[j]:
				(ZTempA[j-1],ZTempA[j]) = (ZTempA[j],ZTempA[j-1])
				(zorderA[j-1],zorderA[j]) = (zorderA[j],zorderA[j-1])
				(LabelsTemp[j-1],LabelsTemp[j]) = (LabelsTemp[j],LabelsTemp[j-1])

	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
	plt.xlabel('Time [years]')
	plt.ylabel('Relative RA [mas]')
	plt.xlim(None, needed_param.epoch[len(needed_param.epoch)-1]+plotMarginX)
    	plt.gca().set_color_cycle(KinematicsWindow.colorsFinal)
    	plt.legend(zorderA,LabelsTemp,loc='upper right',numpoints=1,frameon=True,handletextpad=0,markerscale=1,fontsize = 12)
	plt.title(str(needed_param.source_name), fontsize=15)
	plt.savefig('Plots/RA_with_time.png', bbox_inches='tight')
	plt.savefig('Plots/RA_with_time.ps', bbox_inches='tight')
        plt.savefig('Plots/RA_with_time.pdf', bbox_inches = 'tight')

	zorderA = []
	plt.figure(11)
	for i in xrange(len(components_index)):
		zorderA.append(plt.errorbar(needed_param.epoch,self.yComp[i],yerr=self.yerrComp[i],color=KinematicsWindow.colorsFinal[i],fmt='o'))
	for i in xrange(len(self.rCompNiA)):
		zorderA.append(plt.errorbar(KinematicsWindow.WINcompAutomaticIndexNotInAllepoch[i],self.yCompNiA[i],yerr=self.yerrCompNiA[i],color=KinematicsWindow.colorsFinal[last_comp+i+1],fmt='o'))
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
	#plt.ylim(ymin,ymax)
	plt.xlabel('Time [years]')
	plt.ylabel('Relative DEC [mas]')
	plt.xlim(None, needed_param.epoch[len(needed_param.epoch)-1]+plotMarginX)
	ZTempA = np.asarray(KinematicsWindow.zOrderPlot).copy()
	for i in range(len(ZTempA)):
		for j in range(1,len(ZTemp)-i):
			if ZTempA[j-1] > ZTempA[j]:
				(ZTempA[j-1],ZTempA[j]) = (ZTempA[j],ZTempA[j-1])
				(zorderA[j-1],zorderA[j]) = (zorderA[j],zorderA[j-1])
    	plt.gca().set_color_cycle(KinematicsWindow.colorsFinal)
    	plt.legend(zorderA,LabelsTemp,loc='upper right',numpoints=1,frameon=True,handletextpad=0,markerscale=1,fontsize = 12)
	plt.title(str(needed_param.source_name), fontsize=15)
	plt.savefig('Plots/DEC_with_time.png', bbox_inches='tight')
	plt.savefig('Plots/DEC_with_time.ps', bbox_inches='tight')
	plt.savefig('Plots/DEC_with_time.pdf', bbox_inches = 'tight')

	zorderA = []
	plt.figure(12)
	for i in xrange(len(components_index)):
		zorderA.append(plt.errorbar(needed_param.epoch,self.rComp[i],yerr=self.errrComp[i],color=KinematicsWindow.colorsFinal[i],fmt='o'))
	for i in xrange(len(self.rCompNiA)):
		zorderA.append(plt.errorbar(KinematicsWindow.WINcompAutomaticIndexNotInAllepoch[i],self.rCompNiA[i],yerr=self.errrCompNiA[i],color=KinematicsWindow.colorsFinal[last_comp+i+1],fmt='o'))
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
	#ax.set_ylim(bottom=ymin)
	plt.xlabel('Time [years]')
	plt.ylabel('r [mas]')
	plt.xlim(None, needed_param.epoch[len(needed_param.epoch)-1]+plotMarginX)
	ZTempA = np.asarray(KinematicsWindow.zOrderPlot).copy()
	for i in range(len(ZTempA)):
		for j in range(1,len(ZTemp)-i):
			if ZTempA[j-1] > ZTempA[j]:
				(ZTempA[j-1],ZTempA[j]) = (ZTempA[j],ZTempA[j-1])
				(zorderA[j-1],zorderA[j]) = (zorderA[j],zorderA[j-1])
    	plt.gca().set_color_cycle(KinematicsWindow.colorsFinal)
    	plt.legend(zorderA,LabelsTemp,loc='upper right',numpoints=1,frameon=True,handletextpad=0,markerscale=1,fontsize = 12)
	plt.title(str(needed_param.source_name), fontsize=15)
	plt.savefig('Plots/distance_with_time.png', bbox_inches='tight')
	plt.savefig('Plots/distance_with_time.ps', bbox_inches='tight')
	plt.savefig('Plots/distance_with_time.pdf', bbox_inches = 'tight')

	zorderA = []
	plt.figure(13)
	for i in xrange(len(components_index)):
		zorderA.append(plt.errorbar(needed_param.epoch,self.sizeComp[i],yerr=self.sizeerrComp[i],color=KinematicsWindow.colorsFinal[i],fmt='o'))
	for i in xrange(len(self.rCompNiA)):
		zorderA.append(plt.errorbar(KinematicsWindow.WINcompAutomaticIndexNotInAllepoch[i],self.sizeCompNiA[i],yerr=self.sizeerrCompNiA[i],color=KinematicsWindow.colorsFinal[last_comp+i+1]))
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
	#plt.ylim(ymin,ymax)
	plt.xlabel('Time [years]')
	plt.ylabel('Component Size [mas]')
	plt.xlim(None, needed_param.epoch[len(needed_param.epoch)-1]+plotMarginX)
	ZTempA = np.asarray(KinematicsWindow.zOrderPlot).copy()
	for i in range(len(ZTempA)):
		for j in range(1,len(ZTemp)-i):
			if ZTempA[j-1] > ZTempA[j]:
				(ZTempA[j-1],ZTempA[j]) = (ZTempA[j],ZTempA[j-1])
				(zorderA[j-1],zorderA[j]) = (zorderA[j],zorderA[j-1])
    	plt.gca().set_color_cycle(KinematicsWindow.colorsFinal)
    	plt.legend(zorderA,LabelsTemp,loc='upper right',numpoints=1,frameon=True,handletextpad=0,markerscale=1,fontsize = 12)
	plt.title(str(needed_param.source_name), fontsize=15)
	plt.savefig('Plots/componentSize_with_time.png', bbox_inches='tight')
	plt.savefig('Plots/componentSize_with_time.ps', bbox_inches='tight')
	plt.savefig('Plots/componentSize_with_time.pdf', bbox_inches = 'tight')

	zorderA = []
	plt.figure(14)
	for i in xrange(len(components_index)):
		zorderA.append(plt.errorbar(needed_param.epoch,self.fluxComp[i],yerr=self.fluxerrComp[i],color=KinematicsWindow.colorsFinal[i],fmt='o',ls='--'))
	for i in xrange(len(self.rCompNiA)):
		zorderA.append(plt.errorbar(KinematicsWindow.WINcompAutomaticIndexNotInAllepoch[i],self.fluxCompNiA[i],yerr=self.fluxerrCompNiA[i],color=KinematicsWindow.colorsFinal[last_comp+i+1],fmt='o',ls='--'))
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
        plt.yscale('log')
	plt.xlabel('Time [years]')
	plt.ylabel('Flux [Jy]')
	plt.xlim(None, needed_param.epoch[len(needed_param.epoch)-1]+plotMarginX)
	ZTempA = np.asarray(KinematicsWindow.zOrderPlot).copy()
	for i in range(len(ZTempA)):
		for j in range(1,len(ZTemp)-i):
			if ZTempA[j-1] > ZTempA[j]:
				(ZTempA[j-1],ZTempA[j]) = (ZTempA[j],ZTempA[j-1])
				(zorderA[j-1],zorderA[j]) = (zorderA[j],zorderA[j-1])
    	plt.gca().set_color_cycle(KinematicsWindow.colorsFinal)
    	plt.legend(zorderA,LabelsTemp,loc='upper right',numpoints=1,frameon=True,handletextpad=0,markerscale=1,fontsize = 12)
	plt.title(str(needed_param.source_name), fontsize=15)
	plt.savefig('Plots/componentFlux_with_time.png', bbox_inches='tight')
	plt.savefig('Plots/componentFlux_with_time.ps', bbox_inches='tight')
	plt.savefig('Plots/componentFlux_with_time.pdf', bbox_inches = 'tight')


	plt.figure(9)	
	plt.axis('scaled')
	plt.xlim(self.limplot_x1,self.limplot_x2)
	plt.ylim(self.limplot_y2,self.limplot_y1)
    	plot_components(self.pts_arr,self.x_el_arr,self.x_elH_arr,self.y_elH_arr,self.y_el_arr,needed_param.arr,self.orientation)


	if self.load == False:
		components_index = KinematicsWindow.WINcompIndex
		for i in xrange(0,len(components_index)):
			plot_related_components(self.pts_arr,self.x_el_arr,self.x_elH_arr,self.y_elH_arr,self.y_el_arr,needed_param.arr,components_index[i],KinematicsWindow.colorsFinal[i],self.orientation)
			last_comp = i
			if i==0:
				x = linspace(self.xComp[0][0],self.xComp[0][len(self.xComp[0])-1],100,endpoint=False)
				y = linspace(needed_param.arr[0],needed_param.arr[len(needed_param.arr)-1],100,endpoint=False)
				plt.plot(x,y,color=KinematicsWindow.colorsFinal[0],linewidth=3,ls='solid')
			if self.orientation == 'h':
				if i>0:
					x = linspace(self.xComp[i][0],self.xComp[i][len(self.xComp[i])-1],100,endpoint=False)
					#x = linspace(self.xComp[i][len(self.xComp[i])-1],self.xComp[i][0],100,endpoint=False)
					if np.all(np.asarray(self.xComp[i]) > 0): #add for other jet orientation when needed y = y1 +(x-x1)*(y2-y1)/(x2-x1)
						y = -needed_param.arr[0]+ self.yComp[i][0]- (x-self.xComp[i][0])*(needed_param.arr[len(needed_param.arr)-1]- self.yComp[i][len(self.yComp[i])-1]-needed_param.arr[0]+ self.yComp[i][0])/(self.xComp[i][len(self.xComp[i])-1]-self.xComp[i][0])
					else:
						y = -needed_param.arr[0]+self.yComp[i][0] + (x-self.xComp[i][0])*(-needed_param.arr[len(needed_param.arr)-1]+ self.yComp[i][len(self.yComp[i])-1]+needed_param.arr[0]- self.yComp[i][0])/(self.xComp[i][len(self.xComp[i])-1]-self.xComp[i][0])

					plt.plot(x,y,color=KinematicsWindow.colorsFinal[i],linewidth=3,ls='solid')
			lastCompColor = i
		plot_NiA_components(KinematicsWindow.WINpts_arrNiA,KinematicsWindow.WINx_el_arrNiA,KinematicsWindow.WINy_elH_arrNiA,
				KinematicsWindow.WINx_elH_arrNiA,KinematicsWindow.WINy_el_arrNiA,KinematicsWindow.WINarrNiA,
				KinematicsWindow.WINcompAutomaticIndexNotInAll,KinematicsWindow.colorsFinal,last_comp,
				KinematicsWindow.WINpts_arrNiAtemp,KinematicsWindow.WINx_el_arrNiAtemp,KinematicsWindow.WINy_elH_arrNiAtemp,
				KinematicsWindow.WINx_elH_arrNiAtemp,KinematicsWindow.WINy_el_arrNiAtemp,KinematicsWindow.WINarrNiAtemp,
				KinematicsWindow.WINcompAutomaticIndexNiAnewComp,self.orientation)
		for i in xrange(0,len(self.xCompNiA)):
			if self.orientation == 'h':
				x = linspace(self.xCompNiA[i][0],self.xCompNiA[i][len(self.xCompNiA[i])-1],100,endpoint=False)
				#x = linspace(self.xComp[i][len(self.xComp[i])-1],self.xComp[i][0],100,endpoint=False)
				if np.all(np.asarray(self.xCompNiA[i]) > 0): #add for other jet orientation when needed
					y = -KinematicsWindow.WINarrNiA[i][0]+self.yCompNiA[i][0]- (x-self.xCompNiA[i][0])*(KinematicsWindow.WINarrNiA[i][len(KinematicsWindow.WINarrNiA[i])-1]-self.yCompNiA[i][len(self.yCompNiA[i])-1]-KinematicsWindow.WINarrNiA[i][0]+self.yCompNiA[i][0])/(self.xCompNiA[i][len(self.xCompNiA[i])-1]-self.xCompNiA[i][0])
				else:
					y = -KinematicsWindow.WINarrNiA[i][0]+self.yCompNiA[i][0] + (x-self.xCompNiA[i][0])*(-KinematicsWindow.WINarrNiA[i][len(KinematicsWindow.WINarrNiA[i])-1]+self.yCompNiA[i][len(self.yCompNiA[i])-1]+KinematicsWindow.WINarrNiA[i][0]-self.yCompNiA[i][0])/(self.xCompNiA[i][len(self.xCompNiA[i])-1]-self.xCompNiA[i][0])

				plt.plot(x,y,color=KinematicsWindow.colorsFinal[lastCompColor+i+1],linewidth=3,ls='solid')
		plot_maps(needed_param.realDAT,needed_param.ext,needed_param.arr,needed_param.first_contour,self.orientation)
		plt.title(str(needed_param.source_name), fontsize=35)
		if self.orientation == 'h': 
			plt.xlabel('Relative RA [mas]', fontsize=25)
			ax = plt.gca()
			ax.set_yticks(needed_param.arr)
			ax.set_yticklabels(needed_param.date[::-1])
			ax.tick_params(axis='both',size=20,width=1,which='major', labelsize=25)
		if self.orientation == 'v': 
			plt.ylabel('Relative DEC [mas]', fontsize=25)
			ax = plt.gca()
			ax.set_xticks(needed_param.arr)
			ax.set_xticklabels(needed_param.date[::-1])
			ax.tick_params(axis='both',size=20,width=1,which='major', labelsize=25)
		ax = plt.gcf() # get current figure
		ax.set_size_inches(20, 20)
		plt.savefig('Plots/Identified_components_in_map.png', bbox_inches='tight')#,orientation='landscape')
		plt.savefig('Plots/Identified_components_in_map.ps', bbox_inches='tight')#,orientation='landscape')
		plt.savefig('Plots/Identified_components_in_map.pdf', bbox_inches = 'tight')#,orientation='landscape')

#	plt.figure(15)
#	for i in xrange(1,2):
#		plt.plot(self.rComp[i],-self.psiComp[i]*180/np.pi,linestyle = 'None', marker = 'o',color=self.colors[i])#,xerr =np.asarray(self.errrComp[i]), yerr=np.asarray(self.psierrComp[i]),color=self.colors[i],fmt='o')
#	for i in xrange(2,len(self.rComp)):
#		plt.plot(self.rComp[i],self.psiComp[i]*180/np.pi,linestyle = 'None', marker = 'o',color=self.colors[i])#,xerr =np.asarray(self.errrComp[i]), yerr=np.asarray(self.psierrComp[i]),color=self.colors[i],fmt='o')
#	for i in xrange(len(self.rCompNiA)):
#		plt.plot(self.rCompNiA[i],self.psiCompNiA[i]*180/np.pi,linestyle = 'None', marker = 'o',color=self.colors[last_comp+i+1])#,xerr =np.asarray(self.errrComp[i]),yerr=np.asarray(self.psierrCompNiA[i]),color=self.colors[last_comp+i+1],fmt='o')
	plt.figure(15)
	for i in xrange(1,2):
                print self.psierrComp[i]*180./np.pi
		plt.errorbar(self.rComp[i],-self.psiComp[i]*180/np.pi,xerr=None,yerr = np.asarray(self.psierrComp[i]*180./np.pi)[0], linestyle = 'None', marker = 'o',color=self.colors[i])
	for i in xrange(2,len(self.rComp)):
		plt.errorbar(self.rComp[i],self.psiComp[i]*180/np.pi,xerr=None,yerr = np.asarray(self.psierrComp[i]*180./np.pi)[0],linestyle = 'None', marker = 'o',color=self.colors[i])
	for i in xrange(len(self.rCompNiA)):
		plt.errorbar(self.rCompNiA[i],self.psiCompNiA[i]*180/np.pi,xerr=None,yerr = np.asarray(self.psierrCompNiA[i]*180./np.pi)[0],linestyle = 'None', marker = 'o',color=self.colors[last_comp+i+1])

	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
	#plt.ylim(ymin,ymax)
	plt.xlabel('Distance [mas]')
	plt.ylabel('P.A. [degrees]')
	plt.title(str(needed_param.source_name), fontsize=15)
	plt.savefig('Plots/psi_with_distance.png')
	plt.savefig('Plots/psi_with_distance.pdf')

	#plt.show()
	plt.close('all')

    def Fitting(self):
	self.LineValuesX, self.LineValuesY = [], []
	self.ParabolicValuesX,self.ParabolicValuesY = [], []
	self.FitValuesX, self.FitValuesY = [], []

	self.errLineValuesX, self.errParabolicValuesX = [], []
	self.errLineValuesY, self.errParabolicValuesY = [], []

	self.LineValuesXnia, self.LineValuesYnia = [], []
	self.ParabolicValuesXnia,self.ParabolicValuesYnia = [], []

	self.errLineValuesXnia, self.errParabolicValuesXnia = [], []
	self.errLineValuesYnia, self.errParabolicValuesYnia = [], []

	for i in xrange(1,len(self.xComp)):
		chi2line = probfit.Chi2Regression(line, needed_param.epoch, np.asarray(self.xComp[i]), error=np.asarray(self.xerrComp[i]), weights=None)
		chi2parabolic = probfit.Chi2Regression(parabolic, needed_param.epoch, np.asarray(self.xComp[i]), error=np.asarray(self.xerrComp[i]), weights=None) #np.asarray(xerrComp[i])

		try:
			if np.all(np.asarray(self.xComp[i])<0):
				#write the iminuit line (needed for the fitting)
				mLine=iminuit.Minuit(chi2line,c=needed_param.epoch[0],limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5),m = 0.05,limit_m=(-10.,0.))	
			else:
				mLine=iminuit.Minuit(chi2line,c=needed_param.epoch[0],limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5),m = 0.05,limit_m=(0.,10.))	

			#pfinal, covx = curve_fit(line, needed_param.epoch, np.asarray(self.xComp[i]), [0.5,needed_param.epoch[0]],maxfev=5000)
			mParabolic=iminuit.Minuit(chi2parabolic,c=needed_param.epoch[0],limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5),m = 0.05,limit_m=(-10.,10.), acc=-0.5, tmid=needed_param.tMid, fix_tmid = True)	

			mLine.migrad()
			mLine.hesse()
			print mLine.values  
			self.LineValuesX.append(mLine.values)
			self.errLineValuesX.append(mLine.errors)
			#self.LineValuesX.append(pfinal)

			mParabolic.migrad()
			mParabolic.hesse()
			print mParabolic.values
			self.ParabolicValuesX.append(mParabolic.values)
			self.errParabolicValuesX.append(mParabolic.errors)

			chi2lineValueX = np.sum(((line(needed_param.epoch, mLine.values.get('m'), mLine.values.get('c'))-np.asarray(self.xComp[i]))/np.asarray(self.xerrComp[i]))**2)

			chi2ParabolicValueX = np.sum(((parabolic(needed_param.epoch, mParabolic.values.get('m'), mParabolic.values.get('c'), mParabolic.values.get('acc'),mParabolic.values.get('tmid'))-np.asarray(self.xComp[i]))/np.asarray(self.xerrComp[i]))**2)

			if chi2ParabolicValueX < chi2lineValueX:
				self.FitValuesX.append(mParabolic.values)
				self.errFitValuesX.append(mParabolic.errors)
			else:
				self.FitValuesX.append(mLine.values)
				self.errFitValuesX.append(mLine.errors)

		except RuntimeError:
			print 'Covariance is not valid. May be the last Hesse call failed?'     

	for i in xrange(1,len(self.yComp)):
		chi2line = probfit.Chi2Regression(line, needed_param.epoch, np.asarray(self.yComp[i]), error=np.asarray(self.yerrComp[i]), weights=None)
		chi2parabolic = probfit.Chi2Regression(parabolic, needed_param.epoch, np.asarray(self.yComp[i]), error=np.asarray(self.yerrComp[i]), weights=None) #np.asarray(yerrComp[i])

		try:
			if np.all(np.asarray(self.yComp[i])<0):
			#write the iminuit line (needed for the fitting)
				mLine=iminuit.Minuit(chi2line,c=needed_param.epoch[0],m = 0.05,limit_m=(-10.,0.),limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5))
			else:
				mLine=iminuit.Minuit(chi2line,c=needed_param.epoch[0],m = 0.05,limit_m=(0.,10.),limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5))	
			mParabolic=iminuit.Minuit(chi2parabolic,c=needed_param.epoch[0],limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5),m = 0.05,limit_m=(-10.,10.), acc=0.5, tmid=needed_param.tMid, fix_tmid = True)	

			mLine.migrad()
			mLine.hesse()
			print mLine.values  
			self.LineValuesY.append(mLine.values)
			self.errLineValuesY.append(mLine.errors)

			mParabolic.migrad()
			mParabolic.hesse()
			print mParabolic.values
			self.ParabolicValuesY.append(mParabolic.values)
			self.errParabolicValuesY.append(mParabolic.errors)

			chi2lineValueY = np.sum(((line(needed_param.epoch, mLine.values.get('m'), mLine.values.get('c'))-np.asarray(self.yComp[i]))/np.asarray(self.yerrComp[i]))**2)

			chi2ParabolicValueY = np.sum(((parabolic(needed_param.epoch, mParabolic.values.get('m'), mParabolic.values.get('c'), mParabolic.values.get('acc'),mParabolic.values.get('tmid'))-np.asarray(self.yComp[i]))/np.asarray(self.yerrComp[i]))**2)

			print chi2lineValueY,chi2ParabolicValueY
			if chi2ParabolicValueY < chi2lineValueY:
				self.FitValuesY.append(mParabolic.values)
				self.errFitValuesY.append(mParabolic.errors)
			else:
				self.FitValuesY.append(mLine.values)
				self.errFitValuesY.append(mLine.errors)

		except RuntimeError:
			print 'Covariance is not valid. May be the last Hesse call failed?'     

	#nia plt.
	yNiA = KinematicsWindow.WINcompAutomaticIndexNotInAllepoch 
	#errorbar([i],self.xCompNiA[i],yerr=self.xerrCompNiA[i],color=self.colors[last_comp+i+1],fmt='o')
	for i in xrange(0,len(self.xCompNiA)):
		if len(self.xCompNiA[i]) > 4:
			chi2line = probfit.Chi2Regression(line, yNiA[i], np.asarray(self.xCompNiA[i]), error=np.asarray(self.xerrCompNiA[i]), weights=None)
			chi2parabolic = probfit.Chi2Regression(parabolic, yNiA[i], np.asarray(self.xCompNiA[i]), error=np.asarray(self.xerrCompNiA[i]), weights=None) #np.asarray(xerrComp[i])

			try:
				if np.all(np.asarray(self.xCompNiA[i])<0):
					#write the iminuit line (needed for the fitting)
					mLine=iminuit.Minuit(chi2line,c=needed_param.epoch[0],limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5),m = 0.05,limit_m=(-10.,0.))	
				else:
					mLine=iminuit.Minuit(chi2line,c=needed_param.epoch[0],limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5),m = 0.05,limit_m=(0.,10.))	

				#pfinal, covx = curve_fit(line, needed_param.epoch, np.asarray(self.xComp[i]), [0.5,needed_param.epoch[0]],maxfev=5000)
				mParabolic=iminuit.Minuit(chi2parabolic,c=needed_param.epoch[0],limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5),m = 0.05,limit_m=(-10.,10.), acc=-0.5, tmid=needed_param.tMid, fix_tmid = True)	

				mLine.migrad()
				mLine.hesse()
				print mLine.values  
				self.LineValuesXnia.append(mLine.values)
				self.errLineValuesXnia.append(mLine.errors)
				#self.LineValuesX.append(pfinal)

				mParabolic.migrad()
				mParabolic.hesse()
				print mParabolic.values
				self.ParabolicValuesXnia.append(mParabolic.values)
				self.errParabolicValuesXnia.append(mParabolic.errors)


			except RuntimeError:
				print 'Covariance is not valid. May be the last Hesse call failed?'     

	for i in xrange(0,len(self.xCompNiA)):
		if len(self.xCompNiA[i]) > 4:
			chi2line = probfit.Chi2Regression(line, yNiA[i], np.asarray(self.yCompNiA[i]), error=np.asarray(self.yerrCompNiA[i]), weights=None)
			chi2parabolic = probfit.Chi2Regression(parabolic, yNiA[i], np.asarray(self.yCompNiA[i]), error=np.asarray(self.yerrCompNiA[i]), weights=None) #np.asarray(yerrComp[i])

			try:
				if np.all(np.asarray(self.yCompNiA[i])<0):
				#write the iminuit line (needed for the fitting)
					mLine=iminuit.Minuit(chi2line,c=needed_param.epoch[0],m = 0.05,limit_m=(-10.,0.),limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5))
				else:
					mLine=iminuit.Minuit(chi2line,c=needed_param.epoch[0],m = 0.05,limit_m=(0.,10.),limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5))	
				mParabolic=iminuit.Minuit(chi2parabolic,c=needed_param.epoch[0],limit_c=(needed_param.epoch[0]*0.5,needed_param.epoch[0]*1.5),m = 0.05,limit_m=(-10.,10.), acc=0.5, tmid=needed_param.tMid, fix_tmid = True)	

				mLine.migrad()
				mLine.hesse()
				print mLine.values  
				self.LineValuesYnia.append(mLine.values)
				self.errLineValuesYnia.append(mLine.errors)

				mParabolic.migrad()
				mParabolic.hesse()
				print mParabolic.values
				self.ParabolicValuesYnia.append(mParabolic.values)
				self.errParabolicValuesYnia.append(mParabolic.errors)

			except RuntimeError:
				print 'Covariance is not valid. May be the last Hesse call failed?'     

	epochFit = np.linspace(needed_param.epoch[0],needed_param.epoch[len(needed_param.epoch)-1],100,endpoint=False)

	components_index = self.xComp
	componentsnia = self.yCompNiA
	lastComp = 0 
	plt.figure(2)
	for i in xrange(len(components_index)):		
		plt.errorbar(needed_param.epoch,self.xComp[i],yerr=self.xerrComp[i],color=self.colors[i],fmt='o')
	for i in xrange(len(self.LineValuesX)):		
		#plt.plot(epochFit,line(epochFit,self.LineValuesX[i].get('m'),self.LineValuesX[i].get('c')),color=self.colors[i+1],ls='-')
		plt.plot(epochFit,line(epochFit,self.LineValuesX[i].get('m'),self.LineValuesX[i].get('c')),color=self.colors[i+1],ls='-')
		lastComp = i
	for i in xrange(len(componentsnia)):		
		plt.errorbar(yNiA[i],self.xCompNiA[i],yerr=self.xerrCompNiA[i],color=self.colors[lastComp+i+2],fmt='o')
	for i in xrange(len(self.LineValuesXnia)):		
		plt.plot(epochFit,line(epochFit,self.LineValuesXnia[i].get('m'),self.LineValuesXnia[i].get('c')),color=self.colors[lastComp+i+2],ls='-')
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
	plt.xlabel('Time [years]')
	plt.ylabel('Relative RA [mas]')
	plt.savefig('Plots/RA_with_time.png_fitted_line.png', bbox_inches='tight')
	plt.savefig('Plots/RA_with_time_fitted_line.ps', bbox_inches='tight')
	plt.savefig('Plots/RA_with_time_fitted_line.pdf', bbox_inches = 'tight')

	lastComp = 0
	plt.figure(3)
	for i in xrange(len(components_index)):
		plt.errorbar(needed_param.epoch,self.yComp[i],yerr=self.yerrComp[i],color=self.colors[i],fmt='o')
	for i in xrange(len(self.LineValuesY)):
		plt.plot(epochFit,line(epochFit,self.LineValuesY[i].get('m'),self.LineValuesY[i].get('c')),color=self.colors[i+1],ls='-')
		lastComp = i
	for i in xrange(len(componentsnia)):		
		plt.errorbar(yNiA[i],self.yCompNiA[i],yerr=self.yerrCompNiA[i],color=self.colors[lastComp+i+2],fmt='o')
	for i in xrange(len(self.LineValuesYnia)):		
		plt.plot(epochFit,line(epochFit,self.LineValuesYnia[i].get('m'),self.LineValuesYnia[i].get('c')),color=self.colors[lastComp+i+2],ls='-')
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
	plt.xlabel('Time [years]')
	plt.ylabel('Relative DEC [mas]')
	plt.savefig('Plots/DEC_with_time_fitted_line.png', bbox_inches='tight')
	plt.savefig('Plots/DEC_with_time_fitted_line.ps', bbox_inches='tight')
	plt.savefig('Plots/DEC_with_time_fitted_line.pdf', bbox_inches = 'tight')

	DL = needed_param.DL
	z = needed_param.z
	cte = 4.744*10**8
	c = 3*10**10

	self.muLine = []
	self.betaLine = []
	self.muerrLine = []
	self.betaerrLine = []
	self.tejjLine = []
	self.tejjLineerr = []
	for i in xrange(0,len(self.LineValuesX)):
		muMean = np.sqrt(self.LineValuesX[i].get('m')**2+self.LineValuesY[i].get('m')**2) 
		muMeanerr = np.sqrt(self.errLineValuesX[i].get('m')**2+self.errLineValuesY[i].get('m')**2) 
		self.muLine.append(muMean)
		self.muerrLine.append(muMeanerr)
		betaComp = cte*muMean*needed_param.DL/(c*(1.+needed_param.z))
		betaCompErr = cte*muMeanerr*needed_param.DL/(c*(1.+needed_param.z))
		self.betaLine.append(betaComp)
		self.betaerrLine.append(betaCompErr)
		tejjX = needed_param.tMid - line(epochFit,self.LineValuesX[i].get('m'),self.LineValuesX[i].get('c'))/self.LineValuesX[i].get('m')
		tejjY = needed_param.tMid - line(epochFit,self.LineValuesY[i].get('m'),self.LineValuesY[i].get('c'))/self.LineValuesY[i].get('m')
		tejjAv = (self.LineValuesX[i].get('c')/self.errLineValuesX[i].get('c')**2 + self.LineValuesY[i].get('c')/self.errLineValuesY[i].get('c')**2)/(1./self.errLineValuesX[i].get('c')**2+1./self.errLineValuesY[i].get('c')**2)
		self.tejjLine.append(tejjAv)
		self.tejjLineerr.append((1./self.errLineValuesX[i].get('c')**2+1./self.errLineValuesY[i].get('c')**2)**(-0.5))
	for i in xrange(0,len(self.LineValuesXnia)):
		muMean = np.sqrt(self.LineValuesXnia[i].get('m')**2+self.LineValuesYnia[i].get('m')**2) 
		muMeanerr = np.sqrt(self.errLineValuesXnia[i].get('m')**2+self.errLineValuesYnia[i].get('m')**2) 
		self.muLine.append(muMean)
		self.muerrLine.append(muMeanerr)
		betaComp = cte*muMean*needed_param.DL/(c*(1.+needed_param.z))
		betaCompErr = cte*muMeanerr*needed_param.DL/(c*(1.+needed_param.z))
		self.betaLine.append(betaComp)
		self.betaerrLine.append(betaCompErr)
		tejjX = needed_param.tMid - line(epochFit,self.LineValuesXnia[i].get('m'),self.LineValuesXnia[i].get('c'))/self.LineValuesXnia[i].get('m')
		tejjY = needed_param.tMid - line(epochFit,self.LineValuesYnia[i].get('m'),self.LineValuesYnia[i].get('c'))/self.LineValuesYnia[i].get('m')
		tejjAv = (self.LineValuesXnia[i].get('c')/self.errLineValuesXnia[i].get('c')**2 + self.LineValuesYnia[i].get('c')/self.errLineValuesYnia[i].get('c')**2)/(1./self.errLineValuesXnia[i].get('c')**2+1./self.errLineValuesYnia[i].get('c')**2)
		self.tejjLine.append(tejjAv)
		self.tejjLineerr.append((1./self.errLineValuesXnia[i].get('c')**2+1./self.errLineValuesYnia[i].get('c')**2)**(-0.5))

	print 'mu : ', self.muLine
	print 'mu_err: ', self.muerrLine
	print 'beta: ', self.betaLine
	print 'beta_err: ', self.betaerrLine
	print 'tejj', self.tejjLine
	print 'tejj_err', self.tejjLineerr


	self.mu = []
	self.beta = []
	self.muerr = []
	self.betaerr = []

	for i in xrange(0,len(self.FitValuesX)):
		muMean = np.sqrt(self.FitValuesX[i].get('m')**2+self.FitValuesY[i].get('m')**2) 
		muMeanerr = np.sqrt(self.errFitValuesX[i].get('m')**2+self.errFitValuesY[i].get('m')**2) 
		self.mu.append(muMean)
		self.muerr.append(muMeanerr)
		betaComp = cte*muMean*needed_param.DL/(c*(1.+needed_param.z))
		betaCompErr = cte*muMeanerr*needed_param.DL/(c*(1.+needed_param.z))
		self.beta.append(betaComp)
		self.betaerr.append(betaCompErr)

	#print 'mu : ', self.mu
	#print 'mu_err: ', self.muerr
	#print 'beta: ', self.beta
	#print 'beta_err: ', self.betaerr



	final_txt1 = [self.muLine, self.muerrLine,self.betaLine,self.betaerrLine, self.tejjLine,self.tejjLineerr]
	header = np.array([['#mu (mas/yr) errmu (mas/yr) beta (c) errbeta (c) tejj (yr) errtejj (yr)']])
	final_txt=np.swapaxes(final_txt1, 0,1)
	file_name = 'Velocities.txt'
	saver(file_name, header, final_txt, FORMAT='%1.4f')

	plt.figure(4)
	for i in xrange(len(components_index)):		
		plt.errorbar(needed_param.epoch,self.xComp[i],yerr=self.xerrComp[i],color=self.colors[i],fmt='o')
	for i in xrange(len(self.ParabolicValuesX)):
		plt.plot(epochFit,parabolic(epochFit, self.ParabolicValuesX[i].get('m'), self.ParabolicValuesX[i].get('c'), self.ParabolicValuesX[i].get('acc'),needed_param.tMid),color=self.colors[i+1],ls='-')
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
	plt.xlabel('Time [years]')
	plt.ylabel('Relative RA [mas]')
	plt.savefig('Plots/RA_with_time.png_fitted_parabolic.png', bbox_inches='tight')
	plt.savefig('Plots/RA_with_time_fitted_parabolic.ps', bbox_inches='tight')
	plt.savefig('Plots/RA_with_time_fitted_parabolic.pdf', bbox_inches = 'tight')


	plt.figure(5)
	for i in xrange(len(components_index)):
		plt.errorbar(needed_param.epoch,self.yComp[i],yerr=self.yerrComp[i],color=self.colors[i],fmt='o')
	for i in xrange(len(self.ParabolicValuesY)):
		plt.plot(epochFit,parabolic(epochFit, self.ParabolicValuesY[i].get('m'), self.ParabolicValuesY[i].get('c'), self.ParabolicValuesY[i].get('acc'),needed_param.tMid),color=self.colors[i+1],ls='-')
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
	plt.xlabel('Time [years]')
	plt.ylabel('Relative DEC [mas]')
	plt.savefig('Plots/DEC_with_time.png_fitted_parabolic.png', bbox_inches='tight')
	plt.savefig('Plots/DEC_with_time_fitted_parabolic.ps', bbox_inches='tight')
	plt.savefig('Plots/DEC_with_time_fitted_parabolic.pdf', bbox_inches = 'tight')


	plt.close('all')

	plotMarginX = (np.max(needed_param.epoch)-np.min(needed_param.epoch))*1/3.8
	LabelsTempA = np.asarray(KinematicsWindow.Labels).copy()
	LabelsTemp = LabelsTempA.tolist()
	ZTempA = np.asarray(KinematicsWindow.zOrderPlot).copy()
	zorderA = []


	lastComp = 0
	plt.figure(12)
	for i in xrange(len(components_index)):
		if np.all(np.asarray(self.xComp[i])<0) or np.all(np.asarray(self.yComp[i])<0) :
			yaxis = np.asarray(self.rComp[i]).copy()
		else:
			yaxis =- np.asarray(self.rComp[i]).copy()
		zorderA.append(plt.errorbar(needed_param.epoch,yaxis,yerr=self.errrComp[i],color=self.colors[i],fmt='o'))
		lastComp = i
		if i == 0:
			plt.plot(epochFit,epochFit*0,color=self.colors[i],ls='-')
	for i in xrange(len(self.LineValuesY)):
			tempxt = line(epochFit,self.LineValuesX[i].get('m'),self.LineValuesX[i].get('c'))
			tempyt = line(epochFit,self.LineValuesY[i].get('m'),self.LineValuesY[i].get('c'))
			dt = np.sqrt(tempxt**2+tempyt**2)
			if np.all(np.asarray(self.xComp[i+1])<0) or np.all(np.asarray(self.yComp[i+1])<0) :
				plt.plot(epochFit,dt,color=self.colors[i+1],ls='-')
			else:
				plt.plot(epochFit,-dt,color=self.colors[i+1],ls='-')
	for i in xrange(len(componentsnia)):
		if np.all(np.asarray(self.xCompNiA[i])<0) or np.all(np.asarray(self.yCompNiA[i])<0) :
			yaxis = np.asarray(self.rCompNiA[i]).copy()
		else:
			yaxis =- np.asarray(self.rCompNiA[i]).copy()		
		zorderA.append(plt.errorbar(yNiA[i],yaxis,yerr=self.errrCompNiA[i],color=self.colors[lastComp+i+1],fmt='o'))
	for i in xrange(len(self.LineValuesYnia)):
			tempxt = line(yNiA[i],self.LineValuesXnia[i].get('m'),self.LineValuesXnia[i].get('c'))
			tempyt = line(yNiA[i],self.LineValuesYnia[i].get('m'),self.LineValuesYnia[i].get('c'))
			dt = np.sqrt(tempxt**2+tempyt**2)
			if np.all(np.asarray(self.xCompNiA[i])<0) or np.all(np.asarray(self.yCompNiA[i])<0) :
				plt.plot(yNiA[i],dt,color=self.colors[lastComp+i+1],ls='-')
			else:
				plt.plot(yNiA[i],-dt,color=self.colors[lastComp+i+1],ls='-')
	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
	plt.xlabel('Time [years]')
	plt.ylabel('r [mas]')
	for i in range(len(ZTempA)):
		for j in range(1,len(ZTempA)-i):
			if ZTempA[j-1] > ZTempA[j]:
				(ZTempA[j-1],ZTempA[j]) = (ZTempA[j],ZTempA[j-1])
				(zorderA[j-1],zorderA[j]) = (zorderA[j],zorderA[j-1])
				(LabelsTemp[j-1],LabelsTemp[j]) = (LabelsTemp[j],LabelsTemp[j-1])
	plt.xlim(None, needed_param.epoch[len(needed_param.epoch)-1]+plotMarginX)
    	plt.gca().set_color_cycle(KinematicsWindow.colorsFinal)
    	plt.legend(zorderA,LabelsTemp,loc='upper right',numpoints=1,frameon=True,handletextpad=0,markerscale=1,fontsize = 12)
	plt.title(str(needed_param.source_name), fontsize=15)
	plt.savefig('Plots/distance_with_time_fitted.png', bbox_inches='tight')
	plt.savefig('Plots/distance_with_time_fitted.ps', bbox_inches='tight')
	plt.savefig('Plots/distance_with_time_fitted.pdf', bbox_inches = 'tight')

	plt.show()


    def Load(self):
	#HERE

	self.load = True

	inAll = []
	for filename in sorted(glob.glob('ParametersFiles/Parameters*')):   
		inAll.append(filename)  
	NotinAll = []
	for filename in sorted(glob.glob('ParametersFiles/NiAParameters*')):   
		NotinAll.append(filename)  

	for i in xrange(0,len(inAll)):
		data = np.loadtxt(inAll[i])
		self.rComp.append(data[:,3])
		self.errrComp.append(data[:,4])
		self.psiComp.append(data[:,5])
		self.psierrComp.append(data[:,6])
		self.sizeComp.append(data[:,7])
		self.sizeerrComp.append(data[:,8])
		self.fluxComp.append(data[:,1])
		self.fluxerrComp.append(data[:,2])
		self.tbComp.append(data[:,13])
		self.tberrComp.append(data[:,14])
		self.xComp.append(data[:,9])
		self.xerrComp.append(data[:,10])
		self.yComp.append(data[:,11])
		self.yerrComp.append(data[:,12])

	for i in xrange(0,len(NotinAll)):
		data = np.loadtxt(NotinAll[i])
		KinematicsWindow.WINcompAutomaticIndexNotInAllepoch.append(data[:,0])
		self.rCompNiA.append(data[:,3])
		self.errrCompNiA.append(data[:,4])
		self.psiCompNiA.append(data[:,5])
		self.psierrCompNiA.append(data[:,6])
		self.sizeCompNiA.append(data[:,7])
		self.sizeerrCompNiA.append(data[:,8])
		self.fluxCompNiA.append(data[:,1])
		self.fluxerrCompNiA.append(data[:,2])
		self.tbCompNiA.append(data[:,13])
		self.tberrCompNiA.append(data[:,14])
		self.xCompNiA.append(data[:,9])
		self.xerrCompNiA.append(data[:,10])
		self.yCompNiA.append(data[:,11])
		self.yerrCompNiA.append(data[:,12])


	KinematicsWindow.WINr = np.asarray([self.rComp]).copy()
	KinematicsWindow.WINrNiA = np.asarray([self.rCompNiA]).copy()
	KinematicsWindow.WINrerr = np.asarray([self.errrComp]).copy()
	KinematicsWindow.WINrNiAerr = np.asarray([self.errrCompNiA]).copy()
	KinematicsWindow.WINx = np.asarray([self.xComp]).copy()
	KinematicsWindow.WINxNiA = np.asarray([self.xCompNiA]).copy()
	KinematicsWindow.WINy = np.asarray([self.yComp]).copy()
	KinematicsWindow.WINyNiA = np.asarray([self.yCompNiA]).copy()

	n = len(self.rComp)+len(self.rCompNiA)

	KinematicsWindow.Labels = []

	KinematicsWindow.Labels.append('core')

	for i in xrange(1,n):
		KinematicsWindow.Labels.append('c'+str(len(KinematicsWindow.Labels)))

	KinematicsWindow.colorsFinal = self.colors[0:n]

	KinematicsWindow.zOrderPlot = np.arange(1,len(KinematicsWindow.colorsFinal)+1)



class windowNED():

	"""class ExampleApp(QtGui.QMainWindow, design.Ui_MainWindow):
    def __init__(self, parent=None):
        super(ExampleApp, self).__init__(parent)
        self.setupUi(self)
        self.btnBrowse.clicked.connect(self.doSomething)


    def doSomething(self):
        # Code to replace the main window with a new window
        pass """

	def __init__(self, parent=None):
		#needed_param.source_name = '1258-321'#'PKS1258-321'

		searchInNed = searchNEDnoGUI(needed_param.source_name)
		needed_param.DL = searchInNed[0]
		needed_param.z = searchInNed[1]

		print needed_param.DL

		if needed_param.DL == 0.:
			app = QApplication(sys.argv)
	
			self.w = SearchNEDclass()
			self.w.show()
	
			app.exec_()

			self.w.close()
	
	"""DL = 5489 #mpc #for 2013+370
	z = 0.859
	raDeg = 0.
	decDeg = 0.

	print 'DL =', DL, 'z=',z

	needed_param.DL = DL
	needed_param.z = z"""

class needed_param():

	path = os.getcwd()

	#create the needed directories
	if not os.path.exists('CONV'):
		os.makedirs('CONV')
	if not os.path.exists('SHIFT'):
		os.makedirs('SHIFT')
	if not os.path.exists('Plots'):
		os.makedirs('Plots')
	if not os.path.exists('ParametersFiles'):
		os.makedirs('ParametersFiles')

	#store the uvf, mod and fits files of the maps in a list
	models = []
	for filename in sorted(glob.glob('modelfit/*.mod*')):   
		models.append(filename)  

	cleanmod = []
	for filename in sorted(glob.glob('clean/*.mod*')):   
		cleanmod.append(filename)  

	files = []
	for filename in sorted(glob.glob('clean/*.uvf*')):   
		files.append(filename)  

	fits = []
	for filename in sorted(glob.glob('clean/*.fits*')):   
		fits.append(filename) 

	fitsModelfit = []
	for filename in sorted(glob.glob('modelfit/*.fits*')):   
		fitsModelfit.append(filename) 

	"""
	x1, x2, y1,y2 = edges of the original map in mas
	cent_val_x,cent_val_y = center position of the map in mas
 	rms_maps = rms of the maps
	datamax = peak of the maps
	mapsize = mapsize of each map
	pixel size = cellsize of each map
	"""
	x1, x2, y1,y2, cent_val_x,cent_val_y, rms_maps,datamax = np.array([0.]*len(fits)),np.array([0.]*len(fits)),np.array([0.]*len(fits)),np.array([0.]*len(fits)),np.array([0.]*len(fits)),np.array([0.]*len(fits)),np.array([0.]*len(fits)),np.array([0.]*len(fits))
	mapsize_arr = np.array([0.]*len(fits))
	pixelsize_arr = np.array([0.]*len(fits))
	datamax = np.array([0.]*len(fits))
	rms_fits = np.array([0.]*len(fits))

	#list with the real data of each map 
	realDAT = []

	#order the uvf,mod, fits and modelfit mod with time. The four of them for each epoch need to saved by the user at the same time
	#this ordering is done by looking the files
	ordered = order_by_date(files,cleanmod,models,fits)
	#from the ordered files it gets:
	#the epoch in years for each epoch
	#the circular beam for each epoch
	#and the input list but with the items of the list ordered in time (epoch and beam also ordered)
	epoch = ordered[0]
	bmajReal, bminReal, bpaReal, beam_arr, files,modelsC,fits,modelfit = ordered[2],ordered[3],ordered[4],ordered[5], ordered[6],ordered[7],ordered[8],ordered[9]
	orderedFIT = order_by_date(fitsModelfit,models,models,fitsModelfit)
	modelfit = orderedFIT[7]
	models = modelsC


	#takes the source name from the header
	header = take_header(fits[0])
	source_name = header[8]
	freq = header[5]

	#getting important values of the source as DL and z
	#NedParam = searchNed(source_name)
	
	DL = 0.
	z = 0.

	#takes the date of the observation for each source but in the format existing in the header
	date = []
	for i in xrange(len(fits)):
		header = take_header(fits[i])
		date.append(header[4])


	#calculates the middle time (fist obs + last obs)/2
	tMid = (epoch[0]+epoch[len(epoch)-1])/2.

	#takes the maximum value of the circular beam and select the bmaj, bmin and bpa accordint to it
	beam = np.max(beam_arr)
	bmaj, bmin, bpa = beam, beam, 0

	#from the header of the ordered files it reads the mapsize and the pixelsize of each epoch, it stores it in mapsize_arr and pixelsize_arr
	for i in xrange(0,len(fits)):
		header = take_header(fits[i])
		mapsize_arr[i] = header[10]
		pixelsize_arr[i] = header[0]

	#in the case of slightly different mapsize and cellsize
	#which should be not the case because for this study you want the pixel and mapsize to be as consistent as possible
	#it also should not be the case because the frequency is the same for all the epoch
	#and also the array should be similar
	#so, if the parameters are slightly different, the minimum pixelsize is chosen as the pixelsize for the convolution
	#for the mapsize, the bigger mapsize is chose, it is multiplied by two because difmap stores in the header half of the used value
	mapsize = 2*np.max(mapsize_arr)
	pixelsize = np.min(pixelsize_arr)

	#create the file names for the convolved file for each observation
	conv_files = []
	for i in xrange(0,len(fits)):
		conv_files.append('CONV/'+str('%1.2f' % (epoch[i]))+'convolved_with_beam'+str('%1.2f' % (beam))+'.fits')

	#convolves the files with the previous obtained common values
	for i in xrange(len(files)):
		convolve_difmap([files[i]],[models[i]],bmaj,bmin,bpa,0.,0.,mapsize,pixelsize,[conv_files[i]])
		rms_maps[i] = search_rms()[0]*1000
		os.system('rm difmap.log*\n')

	#reads the modelfit files and obtains the modelfit values and errors

	if os.path.isfile('pos_errors.dat'):
		errors = True
	else:
		errors = None

	mod_parameters = read_modfile(modelfit,beam,bmajReal,bminReal,z,freq,rms_maps/1000.,errors)
	"""
	r, errr = radial distance and error of the component
	psi, errpsi = position angle and error of the component
	size, errsize =  size and error of the component
	"""
	r, errr = mod_parameters[0], mod_parameters[1]     
	psi, errpsi = mod_parameters[2], mod_parameters[3]     
	size, errsize = mod_parameters[4], mod_parameters[5]     
	flux, errflux = mod_parameters[8], mod_parameters[9]    
	tb, errtb = mod_parameters[6], mod_parameters[7]    

	#with the radial distance and position angle, the central positions of the components in RA and DEC are calculated
	"""
	x, errx = position in RA and error of the component
	y, erry = position in DEC and error of the component
	"""
	x_and_y = x_y(r,errr,psi,errpsi,errors)
	x, errx = np.asarray(x_and_y[0]), np.asarray(x_and_y[1])
	y, erry = np.asarray(x_and_y[2]), np.asarray(x_and_y[3])

	#for plotting the components in the map
	"""
	pts_arr = points for drawing the external countour of the ellipse, i.e., the ellipse itself

	x_el_arr = points for the x axis in the x direction of the ellipse. They are between (x_cent_component - size component) and (x_cent_component + size component).They are a total of 50
	y_elH_arr = points for the y axis in the x direction of the ellipse. It is a constant, so it is the same value 50 times, for using it with x_el_arr 

	y_el_arr = points for the y axis in the y direction of the ellipse. They are between (y_cent_component - size component) and (y_cent_component + size component). They are a total of 50
	x_elH_arr = points for the x axis in the y direction of the ellipse. It is a constant, so it is the same value 50 times, for using it with y_el_arr 

	"""
	pts_arr=[]
	pt_arr=[]
	x_el_arr=[]
	x_elH_arr=[]
	y_el_arr=[]
	y_elH_arr=[]

	ellipse_plot = ellipse_axis_lines(x,y,size)
	pts_arr,pt_arr = ellipse_plot[0], ellipse_plot[1]
	x_el_arr,y_el_arr = ellipse_plot[2], ellipse_plot[3]
	x_elH_arr,y_elH_arr = ellipse_plot[4], ellipse_plot[5]  	

	#print x_el_arr
	#print 'second', x_el_arr[0], len(x_el_arr[0]), len(x_el_arr)
	

	#ext = []

	#for the convolved files it reads:
	#the data of each epoch (realDAT)
	#the parameters x1,x2,y1,y2 (edges of the map in mas)
	#the center position of the map in RA and DEC in mas
	#creates array ext, containing the edges coordinates
	#because the mapsize is the same, I only put this one
	for i in xrange(0,len(conv_files)):
		map_data = read_map(fits[i])#conv_files[i])

		realDAT.append(map_data[0])	
		x1[i] = map_data[1]
		x2[i] = map_data[2]
		y1[i] = map_data[3]
		y2[i] = map_data[4]
		cent_val_x[i] = map_data[5]
		cent_val_y[i] = map_data[6]

		ext = np.asarray([x1[i],x2[i],y1[i],y2[i]])

		header = take_header(conv_files[i])

		rms_fits[i] = header[9] 
		datamax[i] = header[11]

	#picks the value of the first contour for plotting the maps
	first_contour = []
	for i in xrange(0,len(rms_maps)):
		first_contour.append(0.12*rms_maps[i])
		#first_contour.append(3*rms_maps[i])
		#first_contour.append(0.2*datamax[i])
		#first_contour.append(200*rms_fits[i])

	#setting the contours		    
	
	#because the epochs are plotted one under the other
	#it calculates the positions to space them evenly in the y axis
	lenght = len(x)
	if len(x)%2==0:
		lims = lenght/2*20
		arr = np.linspace(-lims,lims,lenght,endpoint=False)
	else:
		lims = trunc(lenght/2)*20
		arr = np.linspace(-lims,lims,lenght)


if __name__ == '__main__':
	windowNED()
	app = QApplication(sys.argv)

	w = KinematicsWindow()
	w.show()

	app.exec_()

#windowNED()	
#main()
