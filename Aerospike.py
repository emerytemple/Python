#! /usr/bin/python

from Supersonic import *
import matplotlib.pyplot as plt

def linear_internal(tht, xpe, rrre, gamma, nei, ne):
	''' Linear internal-external aerospike contour
	'''
	rt = 1.0

	xc = []
	yc = []

	xr = []
	yr = []

	me = Supersonic('Mar','x',xpe,'g',gamma)
	ve = Supersonic('nu','M',me,'g',gamma)
	vei = (ve-tht)/2.0

	tht = tht*(pi/180.0)
	vei = vei*(pi/180.0)
	ve = ve*(pi/180.0)

	xc.append(0.0)
	yc.append(0.0)

	xr.append(-rt*sin(tht))
	yr.append(rt*cos(tht))

	thx = tht
	thxc = tht

	dv = (vei-0.0)/(nei-1.0)
	for i in range(1,nei+1):
		vx = (i-1)*dv
		thx = (i-1)*dv

		mx = Supersonic('Mnu','v',vx*(180.0/pi),'g',gamma)
		ax = Supersonic('mu','M',mx)*(pi/180.0)

		angle = 0.5*vx+abs(tht);
		length = sqrt(2*((rrre*rt)**2)*(1.0-cos(vx)))

		xr.append(xr[0]+(length*cos(angle)))
		yr.append(yr[0]+(length*sin(angle)))

		yc.append((((xc[i-1]-xr[i])*tan(thx-ax)*tan(thxc))+(yr[i]*tan(thxc))-(yc[i-1]*tan(thx-ax)))/(tan(thxc)-tan(thx-ax)))
		xc.append((yr[i]-yc[i-1]+(xc[i-1]*tan(thxc))-(xr[i]*tan(thx-ax)))/(tan(thxc)-tan(thx-ax)))

		thxc = thx;

	dv = (ve-vei)/(ne-1.0)
	for i in range(1,ne-1):
		thxc = thxc-dv
		vx = vx + dv

		mx = Supersonic('Mnu','v',vx*(180.0/pi),'g',gamma)
		ax = Supersonic('mu','M',mx)*(pi/180.0)		

		yr.append((((xr[nei+i-1]-xc[nei])*tan(thxc+ax)*tan(thx))+(yc[nei]*tan(thx))-(yr[nei+i-1]*tan(thxc+ax)))/(tan(thx)-tan(thxc+ax)))
		xr.append((yc[nei]-yr[nei+i-1]+(xr[nei+i-1]*tan(thx))-(xc[nei]*tan(thxc+ax)))/(tan(thx)-tan(thxc+ax)))

		thx = thxc

	rexit = yr[len(yr)-1]
	for i in range(1,len(yr)):
		yr[i] = rexit-yr[i]
	for i in range(1,len(yc)):
		yc[i] = rexit-yc[i]

	return xc, yc, xr, yr

def axisymmetric_internal(tht, xpe, rrre, gamma, nei, ne):
	''' Axisymmetric internal-external aerospike contour
	'''

	xc = []
	yc = []

	xr = []
	yr = []
	

	me = Supersonic('Mar','x',xpe,'g',gamma)
	ve = Supersonic('nu','M',me,'g',gamma)

	vei = (ve-tht)/2.0

	mei = Supersonic('Mnu','v',vei,'g',gamma)
	uei = Supersonic('mu','M',mei)*(pi/180.0)
	xpei = Supersonic('A/A*','M',mei,'g',gamma)

	tht = tht*(pi/180.0)
	ve = ve*(pi/180.0)
	vei = vei*(pi/180.0)

	aei = uei+vei-ve

	rpre = sqrt(1.0-mei*xpei*sin(aei)/xpe)
	xpre = (rpre-1.0)*cos(aei)/sin(aei)

	tht = pi/2.0-tht

	dm = (mei-1.0)/(nei-1.0)

	for i in range(1,nei+1):
		mx = 1.0+((i-1)*dm)

		ux = Supersonic('mu','M',mx)*(pi/180.0)
		vx = Supersonic('nu','M',mx,'g',gamma)*(pi/180.0)
		xpx = Supersonic('A/A*','M',mx,'g',gamma)

		bx = tht-(pi/2.0)-vx+abs(vei-ve)
		c = 2.0*rrre*sin(0.5*bx)
		psi = pi-tht+vx-0.5*(pi-bx)

		ax = 2.0*vei-ve-vx+ux

		yr.append(rpre+c*sin(psi))
		xr.append(xpre-c*cos(psi))

		yc.append(sqrt((yr[i-1]**2)+(mx*xpx*sin(ax)/xpe)))
		xc.append(xr[i-1]+(yc[i-1]-yr[i-1])*cos(ax)/sin(ax))

	dm = (me-mei)/(ne-1.0)
	for i in range(1,ne+1):
		mx = mei+((i-1)*dm)

		ux = Supersonic('mu','M',mx)*(pi/180.0)
		vx = Supersonic('nu','M',mx,'g',gamma)*(pi/180.0)
		xpx = Supersonic('A/A*','M',mx,'g',gamma)

		yr.append(sqrt(abs(1.0-mx*xpx*sin(ve-vx+ux)/xpe)))
		xr.append((1.0-yr[len(yr)-1])/tan(ve-vx+ux))

	return xc, yc, xr, yr
'''
#	return {'xc':xc, 'yc':yc, 'xr':xr, 'yr':yr}



	print(' ')
	for i in range(0,len(xr)):
		print(xc[i], yc[i], xr[i], yr[i])

axisymmetric_internal(0.0, 26.0, 9.543908, 1.4, 100, 100)
plt.show()
linear_internal(0.0, 26.0, 9.543908, 1.4, 100, 100)
plt.show()
'''
