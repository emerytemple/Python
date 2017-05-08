#! /usr/bin/python

from Aerospike import *
from Station import *
from Supersonic import *

import matplotlib.pyplot as plt

def main():
	stagnation = Station()
	chamber = Station()
	throat = Station()
	internal = Station()
	external = Station()
	atmosphere = Station()

# inputs

	atmosphere.pressure = 0.032 # Pa
	atmosphere.temperature = 3 # K

	stagnation.pressure = 20*(101325/14.696) # Pa
	stagnation.temperature = ((1000.0-32.0)*(5.0/9.0))+273.15 # K

	chamber.area_ratio = 150.0
	external.area_ratio = 5.0

	total_thrust = 0.125 # N

	gravity = 9.81 # m/s^2
	gamma = 1.4
	R = 287 # J/kg.K

# heat transfer coefficient

D = 0.05
Ar
HD = 2.5
Re = 10000
Pr = 0.712
k = 0.412

K = (1.0+((HD/(0.6/sqrt(Ar)))**6))**(-0.05)
G = 2.0*sqrt(Ar)*((1.0-(2.2*sqrt(Ar)))/(1.0+(0.2*(HD-6)*sqrt(Ar))))
F2 = 0.5*(RE**(2.0/3.0))
Nu = (Pr**0.42)*K*G*F2
h = (Nu*k)/D

# heater
	sigma = 108

	Pow = 2 # W
	V = 3.3 # V

	R = (V*V)/Pow
	I = V/R

	RTF = 1.068
	R68 = R/RTF
	I68 = V/R68

	print(R, I, RTF, R68, I68)

	eps = 0.88
	sbc = 5.67e-8
	H = 0.005
	W = 8
	Acs = H*W
	RL = (sigma/Acs)*0.00001
	L = R68/RL
	P = 2*(H+W)
	SA = P*L
	q = Pow/SA

	print(Acs, RL, L, P, SA, q)

	h = 100
	ts = 255.3722222

	ft = lambda t: (P*h*(t-ts)) + (P*eps*sbc*(t**4-ts**4)) - (I68*I68*RL*1000.0)
	tt = brentq(ft, 1.0, 5000.0)

	conv = P*h*(tt-ts)
	rad = P*eps*sbc*(tt**4-ts**4)
	res = I68*I68*RL*1000.0

	print(conv, rad, res, tt)

# stagnation
	stagnation.mach = 0.0
	stagnation.density = stagnation.pressure/(R*stagnation.temperature)

	stagnation.sonic_velocity = sqrt(gamma*R*stagnation.temperature)
	stagnation.velocity = 0.0

# external
	external.mach = Supersonic('Mar','x',external.area_ratio,'g',gamma)

	tot = Supersonic('To/T','M',external.mach,'g',gamma)
	external.temperature = stagnation.temperature/tot

	pop = Supersonic('Po/P','M',external.mach,'g',gamma)
	external.pressure = stagnation.pressure/pop

	ror = Supersonic('ro/r','M',external.mach,'g',gamma)
	external.density = stagnation.density/ror

	external.sonic_velocity = sqrt(gamma*R*external.temperature)
	external.velocity = external.mach*external.sonic_velocity

	external.nu = Supersonic('nu','M',external.mach,'g',gamma)
	external.mu = Supersonic('mu','M',external.mach)

# internal
	internal.nu = external.nu/2.0
	internal.mach = Supersonic('Mnu','v',internal.nu,'g',gamma)
	internal.mu = Supersonic('mu','M',internal.mach)
	internal.area_ratio = Supersonic('A/A*','M',internal.mach,'g',gamma)

	tot = Supersonic('To/T','M',internal.mach,'g',gamma)
	internal.temperature = stagnation.temperature/tot

	pop = Supersonic('Po/P','M',internal.mach,'g',gamma)
	internal.pressure = stagnation.pressure/pop

	ror = Supersonic('ro/r','M',internal.mach,'g',gamma)
	internal.density = stagnation.density/ror

	internal.sonic_velocity = sqrt(gamma*R*internal.temperature)
	internal.velocity = internal.mach*internal.sonic_velocity

# throat
	throat.mach = 1.0
	throat.area_ratio = 1.0
	
	tot = Supersonic('To/T','M',throat.mach,'g',gamma)
	throat.temperature = stagnation.temperature/tot

	pop = Supersonic('Po/P','M',throat.mach,'g',gamma)
	throat.pressure = stagnation.pressure/pop

	ror = Supersonic('ro/r','M',throat.mach,'g',gamma)
	throat.density = stagnation.density/ror

	throat.sonic_velocity = sqrt(gamma*R*throat.temperature)
	throat.velocity = throat.mach*throat.sonic_velocity

	throat.nu = Supersonic('nu','M',throat.mach,'g',gamma)
	throat.mu = Supersonic('mu','M',throat.mach)

# chamber
	chamber.mach = Supersonic('Marsub','x',chamber.area_ratio,'g',gamma)
	
	tot = Supersonic('To/T','M',chamber.mach,'g',gamma)
	chamber.temperature = stagnation.temperature/tot

	pop = Supersonic('Po/P','M',chamber.mach,'g',gamma)
	chamber.pressure = stagnation.pressure/pop

	ror = Supersonic('ro/r','M',chamber.mach,'g',gamma)
	chamber.density = stagnation.density/ror

	chamber.sonic_velocity = sqrt(gamma*R*chamber.temperature)
	chamber.velocity = chamber.mach*chamber.sonic_velocity

# performance
	a = throat.density*throat.velocity*external.velocity
	b = external.area_ratio*(external.pressure-atmosphere.pressure)
	throat.area = total_thrust/(a+b)

	chamber.area = chamber.area_ratio*throat.area
	internal.area = internal.area_ratio*throat.area
	external.area = external.area_ratio*throat.area

	external.radius = sqrt(external.area/pi)
	print(external.radius)

	mass_flow_rate = throat.density*throat.velocity*throat.area
	momentum_thrust = mass_flow_rate*external.velocity
	pressure_thrust = (external.pressure-atmosphere.pressure)*external.area
	specific_impulse = total_thrust/(gravity*mass_flow_rate)
	characteristic_velocity = (stagnation.pressure*throat.area)/mass_flow_rate # c*
	effective_exhaust_velocity = external.velocity+(((external.pressure-atmosphere.pressure)*external.area)/mass_flow_rate) # c
	thrust_coefficient = total_thrust/(throat.area*stagnation.pressure) # CF

# aerospike contour
	xc, yc, xr, yr = linear_internal(0.0,external.area_ratio,8.0,gamma,100,100)

	for i in range(0,len(xr)):
		if i < len(xc):
			xr[i] = xr[i]*external.radius
			yr[i] = yr[i]*external.radius
			xc[i] = xc[i]*external.radius
			yc[i] = yc[i]*external.radius
		else:
			xr[i] = xr[i]*external.radius
			yr[i] = yr[i]*external.radius

	throat_gap = yc[0]-yr[0]
	print(throat_gap)

	plt.plot(xc, yc, xr, yr)
	plt.show()

# output
	print('')
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10}'.format('','stagnation','chamber','throat','internal','external','atmosphere'))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10}'.format('M', stagnation.mach, chamber.mach, throat.mach, internal.mach, external.mach))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10}'.format('A/A*', '', chamber.area_ratio, throat.area_ratio, internal.area_ratio, external.area_ratio))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:5}'.format('P', stagnation.pressure, chamber.pressure, throat.pressure, internal.pressure, external.pressure,atmosphere.pressure,'Pa'))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:5}'.format('', stagnation.pressure*(14.696/101325.0), chamber.pressure*(14.696/101325.0), throat.pressure*(14.696/101325.0), internal.pressure*(14.696/101325.0), external.pressure*(14.696/101325.0),atmosphere.pressure*(14.696/101325.0),'psi'))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:5}'.format('T', stagnation.temperature, chamber.temperature, throat.temperature, internal.temperature, external.temperature,atmosphere.temperature,'K'))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:5}'.format('rho', stagnation.density, chamber.density, throat.density, internal.density, external.density,'','kg/m^3'))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:5}'.format('a', stagnation.sonic_velocity, chamber.sonic_velocity, throat.sonic_velocity, internal.sonic_velocity, external.sonic_velocity,'','m/s'))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:5}'.format('v', stagnation.velocity, chamber.velocity, throat.velocity, internal.velocity, external.velocity,'','m/s'))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:5}'.format('nu', '', '', throat.nu, internal.nu, external.nu,'','deg'))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:5}'.format('mu', '', '', throat.mu, internal.mu, external.mu,'','deg'))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:5}'.format('A', '', chamber.area, throat.area, internal.area, external.area,'','m^2'))
	print('{0:5} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:5}'.format('', '', chamber.area/(0.0254*0.0254), throat.area/(0.0254*0.0254), internal.area/(0.0254*0.0254), external.area/(0.0254*0.0254),'','in^2'))
	print('')
	print('Performance')
	print('mdot = ', mass_flow_rate, ' kg/s')
	print('F mom = ', momentum_thrust, ' N')
	print('F press = ', pressure_thrust, ' N')
	print('F tot = ', total_thrust, ' N')
	print('Isp = ', specific_impulse, ' s')
	print('c = ', effective_exhaust_velocity, ' m/s')
	print('c* = ', characteristic_velocity, 'm/s')
	print('CF = ', thrust_coefficient)
	print('')
	print('throat gap = ', throat_gap, ' m')
	print('             ', throat_gap/0.0254, ' in')

if __name__ == '__main__':
	main()
