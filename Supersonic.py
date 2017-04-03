#!/usr/bin/python

from scipy.optimize import brentq # rootfinding

from math import *

def Supersonic(string, *arg):
	''' Supersonic isentropic tables.
		All angle given and returned are in degrees.
	'''

	# set inputs
	for i in range(0,len(arg),2):
		if arg[i] == 'b': # angle of oblique shock
			beta1 = arg[i+1]
			beta1 *= (pi/180.0)
		elif arg[i] == 'g': # ratio of specific heats
			gamma = arg[i+1]
			gp1 = gamma+1.0
			gm1 = gamma-1.0
		elif arg[i] == 'M': # mach number
			M = arg[i+1]
		elif arg[i] == 'p': # specific heat
			cp = arg[i+1]
		elif arg[i] == 'PR': # pressure ration (Po/P)
			pr = arg[i+1]
		elif arg[i] == 'R': # ideal gas law constant
			R = arg[i+1]
		elif arg[i] == 't': # new flow angle for oblique shock
			theta = arg[i+1]
			theta *= (pi/180.0)
		elif arg[i] == 'v': # prandtl-meyer function
			nu = arg[i+1]
			nu *= (pi/180.0)
		elif arg[i] == 'x': # area ratio
			xp = arg[i+1]
		else:
			pass

	# ideal gas
		# use try, expect ???

	# 1-D
	if string == 'P*/Po':
		retval = (2.0/gp1)**(gamma/gm1)
	elif string == 'Po/P':
		ToT = 1.0+((gm1/2.0)*(M**2.0))
		retval = ToT**(gamma/gm1)
	elif string == 'r*/ro':
		retval = (2.0/gp1)**(1.0/gm1)
	elif string == 'ro/r':
		ToT = 1.0+((gm1/2.0)*(M**2.0))
		retval = ToT**(1.0/gm1)
	elif string == 'T*/To':
		retval = 2.0/gp1
	elif string == 'To/T':
		retval = 1.0+((gm1/2.0)*(M**2.0))

	elif string == 'Mpr':
		fm = lambda m: pr - Supersonic('Po/P','M',m,'g',gamma)
		retval = brentq(fm, 0.01, 50.0)

	# 1-D w/ heat addition
	elif string == 'P/P*h':
		retval = gp1/(1.0+(gamma*(M**2.0)))
	elif string == 'T/T*h':
		PPs = gp1/(1.0+(gamma*(M**2.0)))
		retval = (M**2.0)*(PPs**2.0)
	elif string == 'r/r*h':
		PPs = gp1/(1.0+(gamma*(M**2.0)))
		retval = (1.0/(M**2.0))*(1.0/PPs)
	elif string == 'Po/Po*h':
		PPs = gp1/(1.0+(gamma*(M**2.0)))
		retval = PPs*(((2.0+(gm1*(M**2.0)))/gp1)**(gamma/gm1))
	elif string == 'To/To*h':
		retval = ((gp1*(M**2.0))/((1.0+(gamma*(M**2.0)))**2.0))*(2.0+(gm1*(M**2.0)))

	# 1-D w/ friction
	elif string == 'T/T*f':
		retval = gp1/(2.0+(gm1*(M**2.0)))
	elif string == 'P/P*f':
		TTs = gp1/(2.0+(gm1*(M**2.0)))
		retval = sqrt(TTs)/M
	elif string == 'r/r*f':
		TTs = gp1/(2.0+(gm1*(M**2.0)))
		retval = sqrt(1/TTs)/M
	elif string == 'Po/Po*f':
		TTs = gp1/(2.0+(gm1*(M**2.0)))
		retval = ((1.0/TTs)**(gp1/(2.0*gm1)))/M
	elif string == '4fL*/D':
		TTs = gp1/(2.0+(gm1*(M**2.0)))
		retval = ((1.0-(M**2.0))/(gamma*(M**2.0)))+((gp1/(2.0*gamma))*log(TTs*(M**2.0)))

	# quasi 1-D
	elif string == 'A/A*':
		ToT = 1.0+((gm1/2.0)*(M**2.0))
		retval = sqrt((1.0/(M**2.0))*(((2.0/gp1)*ToT)**(gp1/gm1)))

	elif string == 'Msub': # xp <= 1.0
		gp1o2gm1 = gp1/(2.0*gm1)
		Togp1 = 2.0/gp1
		gm1o2 = gm1/2.0
		fmf = lambda M: ((M*xp)-((Togp1*(1.0+(gm1o2*(M**2.0))))**gp1o2gm1))	
		retval = brentq(fmf,0.0,1.0)
	elif string == 'Msup': # xp >= 1.0
		gp1o2gm1 = gp1/(2.0*gm1)
		Togp1 = 2.0/gp1
		gm1o2 = gm1/2.0
		fmf = lambda M: ((M*xp)-((Togp1*(1.0+(gm1o2*(M**2.0))))**gp1o2gm1))
		retval = brentq(fmf,1.0,50.0)
	elif string == 'Mar':
		if xp == 1.0:
			retval = 1.0
		else:
			fmf = lambda m: xp - Supersonic('A/A*','M',m,'g',gamma)
			retval = brentq(fmf,1.0,50.0)

	# normal shock
	elif string == 'M2': # M >= 1.0
		ToT = 1.0+((gm1/2.0)*(M**2.0))
		retval = sqrt(ToT/((gamma*(M**2.0))-(gm1/2.0)))
	elif string == 'Po2/Po1':
		PoP = 1.0+(((2.0*gamma)/gp1)*((M**2.0)-1.0))
		ror = (gp1*(M**2.0))/(2.0+(gm1*(M**2.0)))
		ToT = PoP/ror
		del_s = (cp*log(ToT))-(R*log(PoP))
		retval = exp(-del_s/R)
	elif string == 'Po2/P1':
		PoP = 1.0+(((2.0*gamma)/gp1)*((M**2.0)-1.0))
		ror = (gp1*(M**2.0))/(2.0+(gm1*(M**2.0)))
		ToT = PoP/ror
		del_s = (cp*log(ToT))-(R*log(PoP))
		Po2Po1 = exp(-del_s/R)
		TT = 1.0+((gm1/2.0)*(M**2.0))
		Po1P1 = TT**(gamma/gm1)
		retval = Po2Po1*Po1P1
	elif string == 'P2/P1':
		retval = 1.0+(((2.0*gamma)/gp1)*((M**2.0)-1.0))
	elif string == 'r2/r1':
		retval = (gp1*(M**2.0))/(2.0+(gm1*(M**2.0)))
	elif string == 'T2/T1':
		PoP = 1.0+(((2.0*gamma)/gp1)*((M**2.0)-1.0))
		ror = (gp1*(M**2.0))/(2.0+(gm1*(M**2.0)))
		retval = PoP/ror

	# oblique shock
	elif string == 'theta':
		tcb = 2.0*(1.0/tan(beta1))
		top = ((M*sin(beta1))**2.0)-1.0
		bot = ((M**2.0)*(gamma+cos(2.0*beta1)))+2.0
		theta = atan(tcb*(top/bot))
		retval = theta*(180.0/pi)
	elif string == 'betas':
		m2 = ((M**2.0)-1.0)**2.0
		m3 = ((M**2.0)-1.0)**3.0
		gm2 = (gm1/2.0)*(M**2.0)
		gp2 = (gp1/2.0)*(M**2.0)
		gp4 = (gp1/4.0)*(M**4.0)
		t2 = tan(theta)**2.0
		lam = sqrt(m2-(3.0*(1.0+gm2)*(1.0+gp2)*t2))
		xi = (m3-(9.0*(1.0+gm2)*(1.0+gm2+gp4)*t2))/(lam**3.0)
		delta = 0.0
		cx = acos(xi)
		part = ((4.0*pi*delta)+cx)/3.0
		cp = cos(part)
		tb = ((M**2.0)-1.0+(2.0*lam*cp))/(3.0*(1.0+gm2)*tan(theta))
		beta1 = atan(tb)
		retval = beta1*(180.0/pi)
	elif string == 'betaw':
		m2 = ((M**2.0)-1.0)**2.0
		m3 = ((M**2.0)-1.0)**3.0
		gm2 = (gm1/2.0)*(M**2.0)
		gp2 = (gp1/2.0)*(M**2.0)
		gp4 = (gp1/4.0)*(M**4.0)
		t2 = tan(theta)**2.0
		lam = sqrt(m2-(3.0*(1.0+gm2)*(1.0+gp2)*t2))
		xi = (m3-(9.0*(1.0+gm2)*(1.0+gm2+gp4)*t2))/(lam**3.0)
		delta = 1.0
		cx = acos(xi)
		part = ((4.0*pi*delta)+cx)/3.0
		cp = cos(part)
		tb = ((M**2.0)-1.0+(2.0*lam*cp))/(3.0*(1.0+gm2)*tan(theta))
		beta1 = atan(tb)
		retval = beta1*(180.0/pi)

	elif string == 'Mob':
		ftbm = lambda m: (tan(theta) - \
			(2.0*(1.0/tan(beta1))* \
			((((m*sin(beta1))**2.0)-1.0)/ \
			(((m**2.0)*(gamma+cos(2.0*beta1)))+2.0))))
		retval = brentq(ftbm,1.0,50.0)

	# prandtl-meyer expansion waves
	elif string == 'mu':
		mu = asin(1.0/M)
		retval = mu*(180.0/pi)
	elif string == 'nu':
		m2 = (M**2.0)-1.0
		nu = (sqrt(gp1/gm1)*atan(sqrt((gm1/gp1)*m2)))-atan(sqrt(m2))
		retval = nu*(180.0/pi)

	elif string == 'Mnu':
		gp1gm1 = gp1/gm1
		gm1gp1 = gm1/gp1
		fpmf = lambda m: ((sqrt(gp1gm1)*atan(sqrt(gm1gp1*((m**2.0)-1.0))))-atan(sqrt((m**2.0)-1.0))-nu)
		retval = brentq(fpmf, 1.0, 50.0)

	else:
		pass

	return retval
