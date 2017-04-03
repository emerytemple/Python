#!/usr/bin/python

from math import *
from scipy.optimize import brentq
from Supersonic import *

import matplotlib.pyplot as plt

class Aerospike():
	def __init__(self):
			
		self.nr = 10 # number of points on ramp
		self.nc = 10 # number of points on cowl

	def print_external(self, x, y):
		''' Print the contour points.
		'''

		print 'x	y'
		for i in range(40):
			print x[i], y[i]
		print ' '

	def approx_linear_external_contour(self, PR, theta, gamma, throat_len):
		''' Approximate method for the design of 
			an external expansion only 
			two-dimensional plug nozzle.

			"Approximate Method for Plug Nozzle Design"
			by Gianfranco Angelino
			AIAA 1964, Vol. 2, No. 10

			throat is from (0,0) to (x[0],y[0])

			PR is design pressure ratio.
			theta is angle (degrees) throat makes with horizontal.
			gamma is ratio of specific heats.
			throat_len is length of throat.
		'''

		x = []
		y = []

		mach_exit = Supersonic('Mpr','PR',PR,'g',gamma)
		dm = (mach_exit-1.0)/(self.nr-1.0)

		for i in range(self.nr):
			M = 1.0+(i*dm)
			mu = Supersonic('mu','M',M)
			nu = Supersonic('nu','M',M,'g',gamma)

			alpha = (mu+theta-nu)*(pi/180.0)
			xsi = M*Supersonic('A/A*','M',M,'g',gamma)

			# redimensionalize
			L = xsi*throat_len
			x.append(L*cos(alpha))
			y.append(-L*sin(alpha))

		return x, y

	def axisymmetric_external_contour(self, PR, theta, gamma, throat_len):
		''' Approximate method for the design of
			an external expansion only
			axisymmetric plug nozzle.

			"Approximate Method for Plug Nozzle Design"
			by Gianfranco Angelino
			AIAA 1964, Vol. 2, No. 10

			Also Lee & Thompson

			throat is from (0,0) to (x[0],y[0])

			PR is design pressure ratio.
			theta is angle (degrees) throat makes with horizontal.
			gamma is ratio of specific heats.
			throat_len is length of throat.
		'''

		x = []
		y = []

		mach_exit = Supersonic('Mpr','PR',PR,'g',gamma)
		dm = (mach_exit-1.0)/(self.nr-1.0)

		epse = Supersonic('A/A*','M',mach_exit,'g',gamma) # Ae/A*

		alpha = (90.0+theta)*(pi/180.0)
		radius_exit = (throat_len*sin(alpha))/(1.0-sqrt(1.0-(sin(alpha)/epse)))

		for i in range(self.nr):
			M = 1.0+(i*dm)
			mu = Supersonic('mu','M',M)
			nu = Supersonic('nu','M',M,'g',gamma)
			eps = Supersonic('A/A*','M',M,'g',gamma) # A/A*

			alpha = (mu+theta-nu)*(pi/180.0)
			part = abs(1.0-((eps*M*sin(alpha))/epse))
			xsi = (1.0-sqrt(part))/sin(alpha)

			# redimensionalize
			L = xsi*radius_exit
			x.append(L*cos(alpha))
			y.append(-L*sin(alpha))

		return x, y

	def linear_external_contour(self, PR, gamma, throat_len):
		''' Calculate non-dimensional axisymmetric nozzle
			for the desired exit Mach number
			using the Shapiro, Anderson, Denton method
		'''

		x = []
		r = []

		mach_exit = Supersonic('Mpr','PR',PR,'g',gamma)
		nu_exit = Supersonic('nu','M',mach_exit,'g',gamma)*(pi/180.0)

		dnu = 0.01

		x.append(-throat_len*cos((pi/2.0)-nu_exit))
		r.append(throat_len*sin((pi/2.0)-nu_exit))

		npts = 1+int(ceil(nu_exit/dnu))
		for i in range(1,npts):
			theta = nu_exit-(i*dnu)
			nu = min(i*dnu, nu_exit)
			M = Supersonic('Mnu','v',nu*(180.0/pi),'g',gamma)

			LineSlope = Supersonic('mu','M',M,'g',gamma)*(pi/180.0) + max(theta, 0.0)

			a = -tan(theta+dnu)
			b = -tan(LineSlope)
			c = r[i-1] - tan(theta+dnu)*x[i-1]

			r.append((c*b)/(b-a))
			x.append(-c/(b-a))

		xShift = abs(x[0])
		rExit = r[len(r)-1]

		for jj in range(len(x)):
			x[jj] += xShift
			r[jj] = rExit - r[jj]

		# xthroat, rthroat = xShift, rExit

		for i in range(len(x)):
			print x[i], r[i]

		return x, r

	def lee_thompson_internal(self):
		'''

			n1 is number of points
			xp is area ratio
			r is ideal gas law constant (1716)
			te is exit temperature (rankine)
			papc is atmospheric pressure / chamber pressure
			g is gravity (32.17405)
			gamma is ratio of specific heats
		'''

		n1 = 50
		n2 = 50
		peipc = 0.127804525463 # 0.0994474
		xp = 6.18369882353 # 12.0
		rrre = 0.3
		pht = pi/2.0 # sqrt(2.0)/2.0 # 0.993223497048 # pi/2.0
		r = 1718.0
		te = 3000.0
		papc = 0.0151245983948 # 0.01
		g = 32.17405
		gamma = 1.4 # 1.2

		xm = [0 for i in range(n1+2)]
		rx1re = [0 for i in range(n1+2)]
		xx1re = [0 for i in range(n1+2)]
		rx2re = [0 for i in range(n1+2)]
		xx2re = [0 for i in range(n1+2)]
		pxpc = [0 for i in range(n1+2)]


		gp1 = gamma+1.0
		gm1 = gamma-1.0
		exp = gp1/(2.0*gm1)
		c1 = 2.0/gp1
		c2 = gm1/gp1
		c3 = gp1/gm1
		c4 = (1.0-gamma)/gamma
		c5 = 2.0/gm1
		c6 = gamma/gm1

		vvv = lambda m: xp - Supersonic('A/A*','M',m,'g',gamma)
		rm = brentq(vvv,1.0,50.0)
		ve = Supersonic('nu','M',rm,'g',gamma)*(pi/180.0)

		rmei = sqrt(c5*((peipc**c4)-1.0))
		vei = Supersonic('nu','M',rmei,'g',gamma)*(pi/180.0)

		thei = vei - ve
		phei = thei + atan((1/rmei)/sqrt(1-((1/rmei)*(1/rmei))))

		rpre = sqrt(1-(c1*(1+(0.5*gm1*rmei*rmei)))*sin(phei)/xp)
		xpre = (rpre-1)*cos(phei)/sin(phei)
		drm = (rmei-1.0)/n1

		# internal portion of the nozzle
		xm[1] = 1
		for k in range(0,n1+1):
			vx = Supersonic('nu','M',xm[k+1],'g',gamma)*(pi/180.0)
			bx = pht-(pi/2)-vx+abs(thei)
			x1re = 2.0*rrre*sin(0.5*bx)
			psi = pi-pht+vx-0.5*(pi-bx)
			rx1re[k+1] = rpre+x1re*sin(psi)
			xx1re[k+1] = xpre-x1re*cos(psi)

			if 0 == k:
				rx2re[1] = sqrt(rx1re[1]*rx1re[1]+sin(pht)/xp)
				xx2re[1] = xx1re[1]+(rx2re[1]-rx1re[1])*cos(pht)/sin(pht) # use atan
			else:
				phx = 2.0*vei-ve-vx+atan(1.0/(xm[k+1]*sqrt(1.0-(1.0/xm[k+1])**2.0)))
				a2 = c1*(1.0+0.5*gm1*xm[k+1]*xm[k+1])
				rx2re[k+1] = sqrt(rx1re[k+1]*rx1re[k+1]+(a2**(0.5*c3))*sin(phx)/xp)
				xx2re[k+1] = xx1re[k+1]+(rx2re[k+1]-rx1re[k+1])*cos(phx)/sin(phx)

			pxpc[k+1] = (1.0+0.5*gm1*xm[k+1]*xm[k+1])**(-gamma/gm1)
			if k < n1:
				xm[k+2] = xm[k+1] + drm

		print 'internal'
		print 'i	xm	rx1re	xx1re	rx2re	xx2re	pxpc'
		for i in range(1,n1+2):
			print str(i) + '	' + \
				str(xm[i]) + '	' + \
				str(rx1re[i]) + '	' + \
				str(xx1re[i]) + '	' + \
				str(rx2re[i]) + '	' + \
				str(xx2re[i]) + '	' + \
				str(pxpc[i]) + '	'

		# external portion of the nozzle
		xm2 = [0 for i in range(n2+3)]
		rxre = [0 for i in range(n2+2)]
		xxre = [0 for i in range(n2+2)]
		pxpc2 = [0 for i in range(n2+2)]
		sumcg = [0 for i in range(n2+2)]
		sumim = [0 for i in range(n2+2)]
		sumva = [0 for i in range(n2+2)]

		xm2[1] = xm[k+1]
		pxpc2[1] = pxpc[k+1]

		ux = atan(1.0/(xm2[1]*sqrt(1.0-(1.0/xm2[1])**2.0)))
		vx = Supersonic('nu','M',xm2[1],'g',gamma)*(pi/180.0)
		rxre[1] = sqrt(1-((c1*(1+0.5*gm1*xm2[1]*xm2[1]))**exp)*sin(ve-vx+ux)/xp)
		xxre[1] = (1-rxre[1])*cos(ve-vx+ux)/sin(ve-vx+ux)

		ch2 = sqrt((0.5*gp1*xm2[1]*xm2[1])/(1.0+0.5*gm1*xm2[1]*xm2[1]))
		sumcg[1] = gamma*(c1**c6)*ch2*cos(thei)+xp*pxpc2[1]*(1-rxre[1]*rxre[1])

		vt = sqrt(gamma*g*r*(te*(1+0.5*gm1*rm*rm))/(0.5*gp1))
		vc = 1+0.5*gm1*xm2[1]*xm2[1]
		vq = sqrt(gamma*r*g*(te*(1+0.5*gm1*rm*rm))/vc)

		a = (1+0.5*gm1*xm2[1]*xm2[1])**(-c6)
		sumim[1] = vq*(cos(thei)+(1-papc*a*sin(phei)/(xm2[1]*xm2[1]))/gamma)/g
		sumva[1] = vq*(cos(thei)+1/gamma)/g

		drm = (rm-xm2[1])/n2
		xm2[2] = xm2[1]+drm
		pro = pxpc2[1]
		rxo = rxre[1]

		for k in range(1,n2+1):
			ux = atan(1/(xm2[k+1]*sqrt(1-(1/xm2[k+1])**2))) #
			vx = Supersonic('nu','M',xm2[k+1],'g',gamma)*(pi/180.0) #

			rxre[k+1] = c1*(1+0.5*gm1*xm2[k+1]*xm2[k+1])
			rxre[k+1] = rxre[k+1]**(gp1*0.5/gm1)
			rxre[k+1] = sqrt(abs(1-rxre[k+1]*sin(ve-vx+ux)/xp)) # abs is temporary

			xxre[k+1] = (1-rxre[k+1])*cos(ve-vx+ux)/sin(ve-vx+ux) #

			pxpc2[k+1] = (1+0.5*gm1*xm2[k+1]*xm2[k+1])**(-c6) #

			sumcg[k+1] = sumcg[k]+0.5*xp*(pro+pxpc2[k+1])*(rxo*rxo-rxre[k+1]*rxre[k+1]) #

			a = c6
			a = (0.5*gp1)**a
			a = a*vt/(gamma*g) #
			b = 0.5*a*xp #
			sumim[k+1] = sumim[k]+b*(pro+pxpc2[k+1]-2*papc)*(rxo*rxo-rxre[k+1]*rxre[k+1])
			sumva[k+1] = sumva[k]+b*(pro+pxpc2[k+1])*(rxo*rxo-rxre[k+1]*rxre[k+1])

			xm2[k+2] = xm2[k+1] + drm #

			pro = pxpc2[k+1] #
			rxo = rxre[k+1] #

			print ux, vx, a, b, pro, rxo

		print ' '
		print 'external'
		print 'i	xm2	rxre	xxre	pxpc2	sumcg	sumim	sumva'
		for i in range(1,n1+2):
			print str(i) + '	' + \
				str(xm2[i]) + '	' + \
				str(rxre[i]) + '	' + \
				str(xxre[i]) + '	' + \
				str(pxpc2[i]) + '	' + \
				str(sumcg[i]) + '	' + \
				str(sumim[i]) + '	' + \
				str(sumva[i]) + '	'

		return xx1re, rx1re, xx2re, rx2re, xxre, rxre

	def bdenton_axisymmetric_internal(self):
		''' Calculate non-dimensional axisymmetric internal-external nozzle
			for the desired exit Mach number
		'''

		Mei = 2.0
		Mexit = 3.4
		Gamma = 1.4
		rThroat = 1.0
		Beta = 3
		DeltaVAeroD = 0.03

		vmax = Supersonic('nu','M',Mexit,'g',Gamma)*(pi/180.0)
		vRegionOne = Supersonic('nu','M',Mei,'g',Gamma)*(pi/180.0)
		ThetaAeroThroat = vmax - (2*vRegionOne)

		print vmax, vRegionOne, ThetaAeroThroat

		xAeroStreamCowl = [0 for i in range(40)]
		rAeroStreamCowl = [0 for i in range(40)]
		ThetaAeroStreamCowl = [0 for i in range(40)]
		vAeroStreamCowl = [0 for i in range(40)]
		MaAeroStreamCowl = [0 for i in range(40)]
		AlphaAeroStreamCowl = [0 for i in range(40)]

		xAeroStream = [0 for i in range(40)]
		rAeroStream = [0 for i in range(40)]
		ThetaAeroStream = [0 for i in range(40)]
		vAeroStream = [0 for i in range(40)]
		MaAeroStream = [0 for i in range(40)]
		AlphaAeroStream = [0 for i in range(40)]

		xAeroStreamCowl[1] = 0.0
		rAeroStreamCowl[1] = 0.0
		ThetaAeroStreamCowl[1] = ThetaAeroThroat
		vAeroStreamCowl[1] = 0.0
		MaAeroStreamCowl[1] = 1.0
		AlphaAeroStreamCowl[1] = Supersonic('mu','M',MaAeroStreamCowl[1],'g',Gamma)*(pi/180.0)

		ThetaAeroStream[1] = ThetaAeroThroat
		vAeroStream[1] = 0.0
		MaAeroStream[1] = 1.0
		AlphaAeroStream[1] = Supersonic('mu','M',MaAeroStream[1],'g',Gamma)*(pi/180.0)

		if ThetaAeroThroat < 0.0:
			xAeroStream[1] = rThroat*sin(ThetaAeroThroat)
			rAeroStream[1] = rThroat*cos(ThetaAeroThroat)
		elif ThetaAeroThroat > 0.0:
			xAeroStream[1] = -rThroat*sin(ThetaAeroThroat)
			rAeroStream[1] = rThroat*cos(ThetaAeroThroat)
		else:
			xAeroStream[1] = 0.0
			rAeroStream[1] = rThroat

		if ThetaAeroThroat < 0.0:
			xCenter = (rThroat+(Beta*rThroat))*sin(ThetaAeroThroat)
			rCenter = (rThroat+(Beta*rThroat))*cos(ThetaAeroThroat)
		elif ThetaAeroThroat > 0.0:
			xCenter = -(rThroat+(Beta*rThroat))*sin(ThetaAeroThroat)
			rCenter = (rThroat+(Beta*rThroat))*cos(ThetaAeroThroat)
		else:
			xCenter = 0.0
			rCenter = rThroat+(Beta*rThroat)

		ii = 1
		jj = 1
		vRegionOneCheck = 1

		print xAeroStreamCowl[1], rAeroStreamCowl[1], ThetaAeroStreamCowl[1], vAeroStreamCowl[1], MaAeroStreamCowl[1], AlphaAeroStreamCowl[1]
		print xAeroStream[1], rAeroStream[1], ThetaAeroStream[1], vAeroStream[1], MaAeroStream[1], AlphaAeroStream[1]
		print xCenter, rCenter, ii, jj, vRegionOneCheck

		while vRegionOneCheck == 1:
			ii += 1
			jj += 1

			ThetaAeroStream[ii] = ThetaAeroStream[ii-1] + DeltaVAeroD
			vAeroStream[ii] = vAeroStream[ii-1] + DeltaVAeroD
			if vAeroStream[ii] > vRegionOne:
				vAeroStream[ii] = vRegionOne
				ThetaAeroStream[ii] = ThetaAeroThroat + vRegionOne

			vRad = vAeroStream[ii]
			fm = lambda m: vRad*(180.0/pi) - Supersonic('nu','M',m,'g',Gamma)
			MaAeroStream[ii] = brentq(fm, 1.0,50.0)

			AlphaAeroStream[ii] = asin(1/MaAeroStream[ii])

			print ThetaAeroStream[ii], vAeroStream[ii], vRad, MaAeroStream[ii], AlphaAeroStream[ii]

			if ThetaAeroThroat < 0.0:
				TriAngle = ((pi/2)-ThetaAeroThroat)-((pi-vAeroStream[ii])/2)
			elif ThetaAeroThroat > 0.0:
				TriAngle = ((pi/2)+ThetaAeroThroat)-((pi-vAeroStream[ii])/2)
			else:
				TriAngle = (pi/2)-((pi-vAeroStream[ii])/2)

			ChordLength = sqrt(2*((Beta*rThroat)**2)*(1-cos(vAeroStream[ii])))
			DeltaR = ChordLength*sin(TriAngle)
			DeltaX = ChordLength*cos(TriAngle)
			xAeroStream[ii] = xAeroStream[1] + DeltaX
			rAeroStream[ii] = rAeroStream[1] + DeltaR

			print TriAngle, ChordLength, DeltaR, DeltaX, xAeroStream[ii], rAeroStream[ii]

			LineSlope = ThetaAeroStream[ii] - AlphaAeroStream[ii]

			ThetaLast = ThetaAeroStreamCowl[jj-1]
			xLast = xAeroStreamCowl[jj-1]
			rLast = rAeroStreamCowl[jj-1]

			print LineSlope, ThetaLast, xLast, rLast

			a = -tan(ThetaLast)
			b = -tan(LineSlope)
			c = rLast - tan(ThetaLast)*xLast
			d = rAeroStream[ii]-tan(LineSlope)*xAeroStream[ii]
			bot = b - a
			solution1 = ((c*b)-(a*d))/bot
			solution2 = (d-c)/bot

			print a, b, c, d

			rAeroStreamCowl[jj] = solution1
			xAeroStreamCowl[jj] = solution2
			ThetaAeroStreamCowl[jj] = ThetaAeroStream[ii]
			vAeroStreamCowl[jj] = vAeroStream[ii]
			MaAeroStreamCowl[jj] = MaAeroStream[ii]
			AlphaAeroStreamCowl[jj] = AlphaAeroStream[ii]

			print xAeroStreamCowl[jj], rAeroStreamCowl[jj], ThetaAeroStreamCowl[jj], vAeroStreamCowl[jj], MaAeroStreamCowl[jj], AlphaAeroStreamCowl[jj]

			if vAeroStream[ii] >= vRegionOne:
				vRegionOneCheck = 0
			else:
				vRegionOneCheck = 1
			print ' '


		MaContinue = 1
		DeltaVAeroTemp = 0

		while MaContinue == 1:
			ii = ii + 1

			ThetaAeroStreamCowl[jj] = ThetaAeroStreamCowl[jj] - DeltaVAeroTemp
			vAeroStreamCowl[jj] = vAeroStreamCowl[jj] + DeltaVAeroTemp

			vRad = vAeroStreamCowl[jj]
			fm = lambda m: vRad*(180.0/pi) - Supersonic('nu','M',m,'g',Gamma)
			MaAeroStreamCowl[jj] = brentq(fm, 1.0,50.0)

			AlphaAeroStreamCowl[jj] = asin(1/MaAeroStreamCowl[jj])

			print jj, ThetaAeroStreamCowl[jj], vAeroStreamCowl[jj], vRad, MaAeroStreamCowl[jj], AlphaAeroStreamCowl[jj]

			LineSlope = ThetaAeroStreamCowl[jj] + AlphaAeroStreamCowl[jj]

			rStart = rAeroStream[ii-1]
			rIntercept = rAeroStreamCowl[jj] - (tan(LineSlope)*xAeroStreamCowl[jj])
			ThetaLast = ThetaAeroStream[ii-1]
			xLast = xAeroStream[ii-1]
			rLast = rAeroStream[ii-1]

			print LineSlope, rStart, rIntercept, ThetaLast, xLast, rLast

			a = -tan(ThetaLast)
			b = -tan(LineSlope)
			c = rLast - tan(ThetaLast)*xLast
			d = rIntercept
			bot = b - a
			solution1 = ((c*b)-(a*d))/bot
			solution2 = (d-c)/bot

			print a, b, c, d

			rAeroStream[ii] = solution1
			xAeroStream[ii] = solution2
			ThetaAeroStream[ii] = ThetaAeroStreamCowl[jj]
			vAeroStream[ii] = vAeroStreamCowl[jj]
			MaAeroStream[ii] = MaAeroStreamCowl[jj]
			AlphaAeroStream[ii] = AlphaAeroStreamCowl[jj]

			print ii, xAeroStream[ii], rAeroStream[ii], ThetaAeroStream[ii], vAeroStream[ii], MaAeroStream[ii], AlphaAeroStream[ii]

			DeltaVAeroTemp = DeltaVAeroD
			print ' '
			if vAeroStreamCowl[jj] >= vmax:
				MaContinue = 0
			else:
				MaContinue = 1

		rExit = rAeroStream[ii]
		for ll in range(1,ii+1):
			rAeroStream[ll] = rExit - rAeroStream[ll]
		for ll in range(1,jj+1):
			rAeroStreamCowl[ll] = rExit - rAeroStreamCowl[ll]

		rCenter = rExit - rCenter

		print rExit, rAeroStream, rAeroStreamCowl, rCenter

		xi = xAeroStreamCowl
		yi = rAeroStreamCowl
		xe = xAeroStream
		ye = rAeroStream

		return xi, yi, xe, ye

	# rao_optimized_nozzle(self):

class Bell():
	def __init__(self):
		pass

	def ideal_linear(self, dmch, xth1, nch, gamma): # nozflo
		''' Compute the contour of a two-dimensional minimum-length nozzle
			for a design exit Mach number.

			dmch is design Mach number.
			xth1 is angle of first characteristic
			nch is number of characteristics
		'''

		nodes = nch+nch*(nch+1)/2

		xkp = [0 for i in range(nodes+1)]
		xkm = [0 for i in range(nodes+1)]
		xth = [0 for i in range(nodes+1)]
		xnu = [0 for i in range(nodes+1)]
		xm = [0 for i in range(nodes+1)]
		xmu = [0 for i in range(nodes+1)]
		x = [0 for i in range(nodes+1)]
		y = [0 for i in range(nodes+1)]
		slpm = [0 for i in range(nodes+1)]
		slpp = [0 for i in range(nodes+1)]
		dydxm = [0 for i in range(nodes+1)]
		dydxp = [0 for i in range(nodes+1)]

		xth[1] = xth1
	
		dnu = Supersonic('nu','M',dmch,'g',gamma)
		thm = 0.5*dnu
		dth = (thm-xth1)/(nch-1.0)

		xnu[1] = xth[1]
		xkm[1] = xth[1]+xnu[1]
		xkp[1] = xth[1]-xnu[1]

		ib = 1
		ip2 = ib+1

		for i in range(ip2,nch+1):
			xkp[i] = xkp[i-1]
			xth[i] = xth[i-1]+dth
			xnu[i] = xth[i]
			ip = i
			xkm[i] = xth[i]+xnu[i]

		for ncch in range(nch-1,0,-1):
			xkp[ip+1] = xkp[ip]
			xkm[ip+1] = xkm[ip]
			xth[ip+1] = xth[ip]
			xnu[ip+1] = xnu[ip]
			ichm = ip2
			ib = ip + 2
			xkm[ib] = xkm[ip2]
			xth[ib] = 0.0
			xkp[ib] = xth[ib]*2.0-xkm[ib]
			xnu[ib]=(xkm[ib]-xkp[ib])*0.5
			if ncch == 1:
				ip = ip + 2
			ip2 = ib+1
			iend = ib+ ncch-1
			for i in range(ip2,iend+1):
				ichm = ichm+1
				xkp[i]=xkp[i-1]
				xkm[i]=xkm[ichm]
				xth[i]=(xkm[i]+xkp[i])*0.5
				xnu[i]=(xkm[i]-xkp[i])*0.5
				ip=i

		xkp[ip+1] = xkp[ip]
		xkm[ip+1] = xkm[ip]
		xth[ip+1] = xth[ip]
		xnu[ip+1] = xnu[ip]

		ib = 1
		ip2=ib+1
		ncch=nch
		rconv = pi/180.0
		for i in range(1,ip+2):
			# pmf, given nu, find M and mu
			cnu = lambda m: xnu[i] - Supersonic('nu','M',m,'g',1.4)
			xm[i] = brentq(cnu, 1.0, 50.0)
			xmu[i] = Supersonic('mu','M',xm[i],'g',1.4)

			dydxm[i]=tan((xth[i]-xmu[i])*rconv)
			dydxp[i]=tan((xth[i]+xmu[i])*rconv)

		x[1] = -1.0/dydxm[1]
		y[1] = 0.0
		for i in range(2,nch+2):
			isp=i
			slpm[i]=dydxm[i]
			if i == nch+1:
				slpm[i]=tan(xth[i]*rconv)
			slpp[i]=(dydxp[i]+dydxp[i-1])*0.5
			x[i], y[i] = self.simeq(x[i-1],y[i-1],0.0,1.0,slpm[i],slpp[i])

		for ncch in range(nch-1,0,-1):
			ichm=ip2
			ib=isp+1
			slpm[ib] = (dydxm[ib]+dydxm[ichm])*0.5
			x[ib] = -y[ichm]/slpm[ib]+x[ichm]
			y[ib] = 0.0
			ip2 = ib + 1
			iend=ib+ncch
			for i in range(ip2,iend+1):
				ichm=ichm+1
				slpp[i]=(dydxp[i]+dydxp[i-1])*0.5
				slpm[i]=(dydxm[i]+dydxm[ichm])*0.5
				if i == iend:
					slpm[i]=(tan(xth[i]*rconv)+tan(xth[ichm]*rconv))*0.5
				x[i], y[i] = self.simeq(x[i-1],y[i-1],x[ichm],y[ichm],slpm[i],slpp[i])
				isp=i

		print 'i	K-	K+	THETA	NU	M	MU	X	Y'
		for i in range(1,nodes+1):
			print str(i) + '	' + \
				str(xkm[i]) + '	' + \
				str(xkp[i]) + '	' + \
				str(xth[i]) + '	' + \
				str(xnu[i]) + '	' + \
				str(xm[i]) + '	' + \
				str(xmu[i]) + '	' + \
				str(x[i]) + '	' + \
				str(y[i]) + '	'

		return x, y

	def simeq(self, x1, y1, x2, y2, s2, s1):
		xp = (y1-y2+x2*s2-x1*s1)/(s2-s1)
		yp = s1*(xp-x1)+y1
		return xp, yp

	def ideal_axisymmetric(self):
		pass

	def bdenton(self):
		''' Calculate the contour ot an anially symmetric supersonic nozzle
			using a combination of Shapiro, Anderson, and Denton methods
		'''

		Position = [0 for i in range(80)]
		ThetaSAD = [0 for i in range(80)]
		xSAD = [0 for i in range(80)]
		rSAD = [0 for i in range(80)]
		vSAD = [0 for i in range(80)]
		AlphaSAD = [0 for i in range(80)]
		MaSAD = [0 for i in range(80)]

		Mei = 2.0
		Mexit = 3.4
		Gamma = 1.2424
		rThroat = 1.0
		Beta = 3
		DeltaVAeroD = 0.1

		DeltaX = DeltaVAeroD
		NumChar = 1
		PointNum = 0
		ii = 1
		jj = 1
		MaContinue = 0
		while MaContinue == 0:
			PointNum = PointNum+NumChar+1
			for ii in range(PointNum-NumChar,PointNum):
				if ii == PointNum - NumChar:
					print PointNum, NumChar, ii, Position
					Position[ii] = 1
				elif ii == PointNum:
					Position[ii] = 3
				else:
					Position[ii] = 2
			for ii in range(PointNum-NumChar,PointNum):
				if Position[ii] == 1:
					ThetaSAD[ii] = NumChar*DeltaVAeroD
					xSAD[ii] = Beta*rThroat*sin(ThetaSAD[ii])
					rSAD[ii] = rThroat+((Beta*rThroat)*(1-cos(ThetaSAD[ii])))
					vSAD[ii] = ThetaSAD[ii]

					vRad = vSAD[ii]
					fm = lambda m: vRad*(180.0/pi) - Supersonic('nu','M',m,'g',Gamma)
					MaSAD[ii] = brentq(fm, 1.0,50.0)

					AlphaSAD[ii] = asin(1/MaSAD[ii])
				elif Position[ii] == 3:
					ThetaSAD[ii] = 0
					rSAD[ii] = 0

					ThetaSADtemp =ThetaSAD[ii-1]
					rSADtemp = rSAD[ii-1]
					xSADtemp = xSAD[ii-1]
					AlphaSADtemp = AlphaSAD[ii-1]
					MaSADtemp = MaSAD[ii-1]
					vSADtemp = vSAD[ii-1]
					DeltaTheta =1.0
					ThetaLast = 100

					while DeltaTheta >= 1e-10:
						xSAD[ii] = ((rSADtemp-(tan(ThetaSADtemp-AlphaSADtemp)*xSADtemp))/(-tan(ThetaSADtemp-AlphaSADtemp)))-rSAD[ii]
						vSAD[ii] = ThetaSADtemp+vSADtemp-ThetaSAD[ii]+((1/(sqrt((MaSADtemp**2)-1)-(1/tan(ThetaSADtemp))))*((rSAD[ii]-rSADtemp)/rSADtemp))

						vRad = vSAD[ii]
						fm = lambda m: vRad*(180.0/pi) - Supersonic('nu','M',m,'g',Gamma)
						MaSAD[ii] = brentq(fm, 1.0,50.0)
						AlphaSAD[ii] = asin(1/MaSAD[ii])

						DeltaTheta = abs(ThetaLast-ThetaSAD[ii])
						if DeltaTheta > 1e-10:
							ThetaLast = ThetaSAD[ii]
							ThetaSADtemp = (ThetaSADtemp+ThetaSAD[ii])/2
							AlphaSADtemp = (AlphaSADtemp+AlphaSAD[ii])/2
							MaSADtemp = (MaSADtemp + MaSAD[ii])/2
							rSADtemp = (rSADtemp+rSAD[ii])/2
							xSADtemp = (xSADtemp+xSAD[ii])/2
							vSADtemp = (vSADtemp+vSAD[ii])/2
					if MaSAD[ii] >= Mexit:
						MaContinue = 1
					else:
						MaContiue = 0
					NumChar = NumChar+1
				else:
					ThetaSADtempRight = ThetaSAD[ii-1]
					rSADtempRight = rSAD[ii-1]
					xSADtempRight = xSAD[ii-1]
					AlphaSADtempRight = AlphaSAD[ii-1]
					MaSADtempRight = MaSAD[ii-1]
					vSADtempRight = vSAD[ii-1]

					ThetaSADtempLeft = ThetaSAD[ii-NumChar]
					rSADtempLeft = rSAD[ii-NumChar]
					xSADtempleft = xSAD[ii-NumChar]
					AlphaSADtempLeft = AlphaSAD[ii-NumChar]
					MaSADtempLeft = MaSAD[ii-NumChar]
					vSADtempleft = vSAD[ii-NumChar]

					DeltaTheta = 1.0
					ThetaLast = 100
					while DeltaTheta >= 1e-10:
						a = tan(ThetaSADtempRight - AlphaSADtempRight)
						b = tan(ThetaSADtempLeft + AlphaSADtempLeft)
						c = rSADtempRight - (a*xSADtempRight)
						d = rSADtempLeft - (b*xSADtempLeft)
						bot = a - b
						solution1 = ((a*d)-(b*c))/bot
						solution2 = (d-c)/bot

						rSAD[ii] = solution1
						xSAD[ii] = solution2

						c = (ThetaSADtempRight+vSADtempRight)+((1/((sqrt(MaSADtempRight**2-1))-(1/tan(ThetaSADtempRight))))*((rSAD[ii]-rSADtempRight)/rSADtempRight))
						if ThetaSADtempLeft == 0.0:
							d = (2*ThetaSADtempLeft)-vSADtempLeft
							bot = 3
						else:
							d = (ThetaSADtempLeft-vSADtempLeft)-((1/((sqrt(MaADtempLeft**2-1))+(1/tan(ThetaSADtempLeft))))*((rSAD[ii]-rSADtempLeft)/rSADtempLeft))
							bot = 2
						solution1 = (c+d)/bot
						solution2 = (c-d)/bot
						ThetaSAD[ii] = solution1
						vSAD[ii] = solution2

						vRad = vSAD[ii]
						fm = lambda m: vRad*(180.0/pi) - Supersonic('nu','M',m,'g',Gamma)
						MaSAD[ii] = brentq(fm, 1.0,50.0)
						AlphaSAD[ii] = asin(1/MaSAD[ii])

						DeltaTheta = abs(ThetaLast-ThetaSAD[ii])
						if DeltaTheta > 1e-10:
							ThetaLast = ThetaSAD[ii]

							ThetaSADtempRight = (ThetaSADtempRight+ThetaSAD[ii])/2
							AlphaSADtempRight = (AlphaSADtempRight+AlphaSAD[ii])/2
							MaSADtempRight = (MaSADtempRight+MaSAD[ii])/2
							rSADtempRight = (rSADtempRight+rSAD[ii])/2
							xSADtempRight = (xSADtempRight+xSAD[ii])/2
							vSADtempRight = (vSADtempRight+vSAD[ii])/2

							ThetaSADtempLeft = (ThetaSADtempLeft+ThetaSAD[ii])/2
							AlphaSADtempLeft = (AlphaSADtempLeft+AlphaSAD[ii])/2
							MaSADtempLeft = (MaSADtempLeft+MaSAD[ii])/2
							rSADtempLeft = (rSADtempLeft+rSAD[ii])/2
							xSADtempLeft = (xSADtempLeft+xSAD[ii])/2
							vSADtempLeft = (vSADtempLeft+vSAD[ii])/2
			
		ii = PointNum+1
		ThetaSAD[ii] = ThetaSAD[ii-(NumChar-1)]
		a = tan(ThetaSAD[ii-(NumChar-1)]+AlphaSAD[ii-(NumChar-1)])
		b = (ThetaSAD[ii-NumChar]+ThetaSAD[ii-(NumChar-1)])/2
		c = rSAD[ii-(NumChar-1)]-(a*xSAD[ii-(NumChar-1)])
		d = rSAD[ii-NumChar]-(b*xSAD[ii-NumChar])
		bot = a - b
		solution1 = ((a*d)-(b*c))/bot
		solution2 = (d-c)/bot

		rSAD[ii] = solution1
		xSAD[ii] = solution2

		for jj in range(ii+1,PointNum-(NumChar-1)):
			ThetaSAD[jj] = ThetaSAD[jj-(NumChar-1)]
			a = tan(ThetaSAD[jj-(NumChar-1)]+AlphaSAD[jj-(NumChar-1)])
			b = (ThetaSAD[jj-NumChar]+ThetaSAD[jj-(NumChar-1)])/2
			c = rSAD[jj-(NumChar-1)]-(a*xSAD[jj-(NumChar-1)])
			d = rSAD[jj-NumChar]-(b*xSAD[jj-NumChar])
			bot = a - b
			solution1 = ((a*d)-(b*c))/bot
			solution2 = (d-c)/bot

			rSAD[jj] = solution1
			xSAD[jj] = solution2

		jj = 1
		Char = 1
		for ii in range(1,NumChar-1):
			xSADcontour[ii] = xSAD[jj]
			rSADcontour[ii] = rSAD[jj]
			jj = jj+(Char+1)
			Char = Char+ 1
		jj = PointNum+1
		for ii in range(NumChar,(NumChar-1)+(NumChar-1)):
			xSADcontour[ii] = xSAD[jj]
			rSADcontour[ii] = rSAD[jj]

class Performance():
	def Calc_TIsp(self,x, y, Po, gamma):
		''' Calculate the Thrust and Isp.
			Throat is at (0,0).
		'''

		Thrust = 0.0
		Isp = 0.0

		length = len(x)

		Pressure = []
		for i in range(length):
			if i == 0:
				Aratio = 1.01
			else:
				Aratio = sqrt((x[i]**2)+(y[i]**2))
			Mach = Supersonic('Mar','x',Aratio,'g',gamma)
			Pressure = Supersonic('Po/P','M',Mach,'g',gamma)
			print Aratio, Mach, Pressure

		return Thrust, Isp

aspike = Aerospike()
x, y = aspike.approx_linear_external_contour(66.1174580574,56.9075144941,1.4, 1.0)
perf = Performance()
T, I = perf.Calc_TIsp(x, y,66.1.4)
#aspike.print_external(x, y)
plt.plot(x, y)
plt.show()

