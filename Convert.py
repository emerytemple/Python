#!/usr/bin/python

from math import pi

def Convert(number, from_unit, to_unit):
	''' converts a number from one measurement system to another
		number - value in from_units to convert
		from_unit - units for the number
		to_unit - units for the result
	'''

	if from_unit in ['g', 'kg', 'sl', 'lbm']: # mass
		number = mass_convert(number, from_unit, to_unit)
	elif from_unit in ['mm', 'cm', 'm', 'km', 'in', 'ft', 'yd', 'mi']: # distance
		number = distance_convert(number, from_unit, to_unit)
	elif from_unit in ['s', 'sec', 'min', 'hr']: # time
		number = time_convert(number, from_unit, to_unit)
	elif from_unit in ['pa', 'atm', 'psi']: # pressure
		number = pressure_convert(number, from_unit, to_unit)
	elif from_unit in ['n', 'lbf']: # force
		number = force_convert(number, from_unit, to_unit)
	elif from_unit in ['c', 'f', 'k', 'r']: # temperature
		number = temperature_convert(number, from_unit, to_unit)
	elif from_unit in ['deg', 'rad']: # angle
		number = angle_convert(number, from_unit, to_unit)
	else:
		pass

	return number

# mass
def mass_convert(number, from_unit, to_unit):
	if from_unit == 'g': # gram
		if to_unit == 'kg':
			number *= 0.001
		elif to_unit == 'sl':
			number *= 6.85218e-5
		elif to_unit == 'lbm':
			number *= 0.00220462
		else:
			pass
	elif from_unit == 'kg': # kilogram
		if to_unit == 'g':
			number *= 1000.0
		elif to_unit == 'sl':
			number *= 0.0685218
		elif to_unit ==	'lbm':
			number *= 2.20462
		else:
			pass
	elif from_unit == 'sl': # slug
		if to_unit == 'g':
			number *= 14593.9
		elif to_unit == 'kg':
			number *= 14.5939
		elif to_unit == 'lbm':
			number *= 32.174
		else:
			pass
	elif from_unit == 'lbm': # pound mass
		if to_unit == 'g':
			number *= 453.592
		elif to_unit == 'kg':
			number *= 0.453592
		elif to_unit == 'sl':
			number *= 0.031081
		else:
			pass

	return number

# distance
def distance_convert(number, from_unit, to_unit):
	if from_unit == 'mm': # millimeter
		if to_unit == 'cm':
			number *= 0.1
		elif to_unit == 'm':
			number *= 0.001
		elif to_unit == 'km':
			number *= 1.0e-6
		elif to_unit == 'in':
			number *= 0.0393701
		elif to_unit == 'ft':
			number *= 0.00328084
		elif to_unit == 'yd':
			number *= 0.00109361
		elif to_unit == 'mi':
			number *= 6.21371e-7
		else:
			pass
	elif from_unit == 'cm': # centimeter
		if to_unit == 'mm':
			number *= 10.0
		elif to_unit == 'm':
			number *= 0.01
		elif to_unit == 'km':
			number *= 1.0e-5
		elif to_unit == 'in':
			number *= 0.393701
		elif to_unit == 'ft':
			number *= 0.0328084
		elif to_unit == 'yd':
			number *= 0.0109361
		elif to_unit == 'mi':
			number *= 6.21371e-6
		else:
			pass
	elif from_unit == 'm': # meter
		if to_unit == 'mm':
			number *= 1000.0
		elif to_unit == 'cm':
			number *= 100.0
		elif to_unit == 'km':
			number *= 0.001
		elif to_unit == 'in':
			number *= 38.3701
		elif to_unit == 'ft':
			number *= 3.28084
		elif to_unit == 'yd':
			number *= 1.09361
		elif to_unit == 'mi':
			number *= 0.000621371
		else:
			pass
	elif from_unit == 'km': # kilometer
		if to_unit == 'mm':
			number *= 1.0e6
		elif to_unit == 'cm':
			number *= 100000.0
		elif to_unit == 'm':
			number *= 1000.0
		elif to_unit == 'in':
			number *= 39370.1
		elif to_unit == 'ft':
			number *= 3280.84
		elif to_unit == 'yd':
			number *= 1093.61
		elif to_unit == 'mi':
			number *= 0.621371
		else:
			pass
	elif from_unit == 'in': # inch
		if to_unit == 'mm':
			number *= 25.4
		elif to_unit == 'cm':
			number *= 2.54
		elif to_unit == 'm':
			number *= 0.0254
		elif to_unit == 'km':
			number *= 2.54e-5
		elif to_unit == 'ft':
			number *= 0.0833333
		elif to_unit == 'yd':
			number *= 0.0277778
		elif to_unit == 'mi':
			number *= 1.5783e-5
		else:
			pass
	elif from_unit == 'ft': # foot
		if to_unit == 'mm':
			number *= 304.8
		elif to_unit == 'cm':
			number *= 30.48
		elif to_unit == 'm':
			number *= 0.3048
		elif to_unit == 'km':
			number *= 0.0003048
		elif to_unit == 'in':
			number *= 12.0
		elif to_unit == 'yd':
			number *= 0.333333
		elif to_unit == 'mi':
			number *= 0.000189394
		else:
			pass
	elif from_unit == 'yd': # yard
		if to_unit == 'mm':
			number *= 914.4
		elif to_unit == 'cm':
			number *= 91.44
		elif to_unit == 'm':
			number *= 0.9144
		elif to_unit == 'km':
			number *= 0.0009144
		elif to_unit == 'in':
			number *= 36.0
		elif to_unit == 'ft':
			number *= 3.0
		elif to_unit == 'mi':
			number *= 0.000568182
		else:
			pass
	elif from_unit == 'mi': # mile
		if to_unit == 'mm':
			number *= 1.60934e6
		elif to_unit == 'cm':
			number *= 160934
		elif to_unit == 'm':
			number *= 1609.34
		elif to_unit == 'km':
			number *= 1.60934
		elif to_unit == 'in':
			number *= 63360.0
		elif to_unit == 'ft':
			number *= 5280.0
		elif to_unit == 'yd':
			number *= 1760.0
		else:
			pass
	
	else:
		pass

	return number

# time
def time_convert(number, from_unit, to_unit):
	if from_unit in ['s', 'sec']: # second
		if to_unit == 'min':
			number *= 0.0166667
		elif to_unit == 'hr':
			number *= 0.000277778
		else:
			pass
	elif from_unit == 'min': # minute
		if to_unit in ['s', 'sec']:
			number *= 60.0
		elif to_unit == 'hr':
			number *= 0.0166667
		else:
			pass
	elif from_unit == 'hr': # hour
		if to_unit in ['s', 'sec']:
			number *= 3600.0
		elif to_unit == 'min':
			number *= 60.0
		else:
			pass
	else:
		pass

	return number

# pressure
def pressure_convert(number, from_unit, to_unit):
	if from_unit == 'pa': # pascal
		if to_unit == 'atm':
			number *= 9.86923e-6
		elif to_unit == 'psi':
			number *= 0.000145038
		else:
			pass
	elif from_unit == 'atm': # atmostphere
		if to_unit == 'pa':
			number *= 101325.0
		elif to_unit == 'psi':
			number *= 14.6959
		else:
			pass
	elif from_unit == 'psi': # pound per square inch
		if to_unit == 'pa':
			number *= 6894.76
		elif to_unit == 'atm':
			number *= 0.068046
		else:
			pass
	else:
		pass

	return number

# force
def force_convert(number, from_unit, to_unit):
	if from_unit == 'n': # newton
		if to_unit == 'lbf':
			number *= 0.224809
		else:
			pass
	elif from_unit == 'lbf': # pound force
		if to_unit == 'n':
			number *= 4.44822
		else:
			pass
	else:
		pass

	return number
	
# temperature
def temperature_convert(number, from_unit, to_unit):
	if from_unit == 'c':
		if to_unit == 'f':
			number = (1.8*number)+32.0
		elif to_unit == 'k':
			number += 273.15
		elif to_unit == 'r':
			number = 1.8*(number+273.15)
		else:
			pass
	elif from_unit == 'f':
		if to_unit == 'c':
			number = (5.0/9.0)*(number-32.0)
		elif to_unit == 'k':
			number = (5.0/9.0)*(number+459.67)
		elif to_unit == 'r':
			number += 459.67
		else:
			pass
	elif from_unit == 'k':
		if to_unit == 'c':
			number -= 273.15
		elif to_unit == 'f':
			number = (1.8*number)-459.67
		elif to_unit == 'r':
			number *= 1.8
		else:
			pass
	elif from_unit == 'r':
		if to_unit == 'c':
			number = (5.0/9.0)*(number-491.67)
		elif to_unit == 'f':
			number -= 459.67
		elif to_unit == 'k':
			number *= 5.0/9.0
		else:
			pass
	else:
		pass

	return number

# angle
def angle_convert(number, from_unit, to_unit):
	if from_unit == 'deg':
		if to_unit == 'rad':
			number *= pi/180.0
		else:
			pass
	elif from_unit == 'rad':
		if to_unit == 'deg':
			number *= 180.0/pi
		else:
			pass
	else:
		pass

	return number

aaa = Convert(1000.0, 'g', 'kg')
bbb = Convert(5280.0, 'ft', 'mi')
ccc = Convert(45.0, 's', 'min')
ddd = Convert(14.6959, 'psi', 'pa')
eee = Convert(0.125, 'n', 'lbf')
fff = Convert(110, 'r', 'k')
ggg = Convert(90.0, 'deg', 'rad')
print(aaa, bbb, ccc, ddd, eee, fff, ggg)

'''

	# energy
	j	joule
	cal	calorie # big or little?
	btu	british thermal unit

	# volume
	ft3	cubic feet
	in3	cubic inch
	m3	cubic meter

	# area
	ft2 square feet
	in2	square inche
	m2	square meter
	mi2	square mile

	# speed
	m/s	meter per second
	mph	mile per hour
'''	
