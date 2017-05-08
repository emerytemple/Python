#!/usr/bin/python

class Airfoil:
	def __init__(self, name):
		if name == '0012':
			pass

class Wing:
	def__init__(self):
		self.cr = 1.0
		self.ct = 1.0
		self.b = 1.0
		self.sweep = 0.0

def main():
	naca0012 = Airfoil('0012')
	wingtest = Wing()
	print(wingtest.cr)

if __name__ == '__main__':
	main()
