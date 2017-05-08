#!/usr/bin/python3.1

import logging
import sys

import ga

import gachrome as gac
import gameiosis as gam
import gaparent as gap
import gaparam as param

def main(args):
	logging.basicConfig(format='', level=logging.DEBUG)

	chrome = gac.Chromosome()
	chrome.pop_chrome(param.num_bool, \
							  param.num_int, \
							  param.num_dbl, \
							  param.num_perm, \
							  param.int_low_bnd, \
							  param.int_high_bnd, \
							  param.dbl_low_bnd, \
							  param.dbl_high_bnd, \
							  param.perm_len)
	chrome.evaluate()
	print(chrome.__str__())

'''
	genalg = ga.GA()
	genalg.evolve()
'''
if __name__ == '__main__':
	main(sys.argv)
