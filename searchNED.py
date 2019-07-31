import urllib2
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

def searchNEDnoGUI(source_name):

	if source_name.find('+') > 0:
		splitedName = source_name.split('+')
		name1 = splitedName[0]
		name2 = splitedName[1]
		response = urllib2.urlopen('https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname='+str(name1)+'%2B'+str(name2)+'&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES')
	if source_name.find('-') > 0:
		splitedName = source_name.split('-')
		name1 = splitedName[0]
		name2 = splitedName[1]
		response = urllib2.urlopen('https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname='+str(name1)+'-'+str(name2)+'&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES')

	html = response.read()

	DLpos = html.find('Luminosity Distance')
	DLpos1st = DLpos+29
	afterDLpos = html[DLpos:]
	DLposEnd = DLpos+afterDLpos.find('Mpc')-2

	#print html[DLpos1st], html[DLposEnd]

	DL= ''
	for i in xrange (DLpos1st,DLposEnd+1):
		DL= DL + html[i]

	try:		
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
		if source_name.find('+') > 0:
			decDeg = decD + decM/(60) + decS/(60*60)
		if source_name.find('-') > 0:
			decDeg = decD - decM/(60) - decS/(60*60)

		#c = SkyCoord(ra=raHours*u.degree,dec=decDays*u.degree)
	
		print 'RA =', ra, 'DEC =', dec
		print 'RA =', raDeg, 'degrees, DEC =', decDeg, 'degrees' 

	except ValueError:
		DL = 0.
		z = 0.

	return DL,z


