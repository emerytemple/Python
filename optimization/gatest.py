#!/usr/bin/python3.1

import logging
import unittest

import gachrome as gac

class ChromosomeTest(unittest.TestCase):

	def test_bitwise_mutation(self):

		# zero length
		testchrome = gac.Chromosome()

		pm = None
		testchrome.bitwise_mutation(pm)

		self.assertEqual(testchrome.bool_chrome, [])

		# change none
		testchrome = gac.Chromosome()

		num_bool = 10
		testchrome.pop_chrome(num_bool,0,0,0,0,0,0,0)
		beforetest = list(testchrome.bool_chrome)

		pm = 0.0
		testchrome.bitwise_mutation(pm)

		self.assertEqual(beforetest, testchrome.bool_chrome)

		# change all
		testchrome = gac.Chromosome()

		num_bool = 15
		testchrome.pop_chrome(num_bool,0,0,0,0,0,0,0)
		beforetest = list(testchrome.bool_chrome)
		ans = [False if i else True for i in beforetest]

		pm = 1.0
		testchrome.bitwise_mutation(pm)

		self.assertEqual(testchrome.bool_chrome, ans)

	def test_creep_mutation(self):

		num_int = 5

		ilow = [-30, -20, -10, 0, 10]
		ihigh = [-10, 0, 10, 20, 30]

		sigma = 1.0

		testchrome = gac.Chromosome()

		# zero length
		pm = 1.0
		testchrome.creep_mutation(pm,[],[],sigma)
		self.assertEqual(testchrome.int_chrome, [])

		# change none
		testchrome.pop_chrome(0,num_int,0,0,ilow,ihigh,0,0)
		before = list(testchrome.int_chrome)

		pm = 0.0
		testchrome.creep_mutation(pm,ilow,ihigh,sigma)
		self.assertEqual(before, testchrome.int_chrome)

		# change all
		testchrome.pop_chrome(0,num_int,0,0,ilow,ihigh,0,0)
		before = list(testchrome.int_chrome)

		pm = 1.0
		testchrome.creep_mutation(pm,ilow,ihigh,sigma)
		for i, chrome in enumerate(testchrome.int_chrome):
			self.assertLessEqual(chrome, ihigh[i])
			self.assertGreaterEqual(chrome, ilow[i])

	def test_insert_mutation(self):

		num_perm = 10

		testchrome = gac.Chromosome()

		# zero length
		testchrome.insert_mutation()
		self.assertEqual(testchrome.perm_chrome, [])

		testchrome.perm_chrome = [0,1,2]
		print('zzzzzzzzzzzzzzzzzzzzzzzz', len(testchrome.permchrome))
		testchrome.insert_mutation()

logging.basicConfig(format='', level=logging.DEBUG)
unittest.main()
