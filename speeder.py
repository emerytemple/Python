#!/usr/bin/python

import math # for sqrt
import random
import sys # handle command line arguments
import unittest

class GA:
	def __init__(self):
		self.pop_size = 10
		self.max_iter = 1

		self.parent_choice = "f"

		self.s = 2.0

		self.num_bool = 5
		self.num_int = 5
		self.num_dbl = 5
		self.num_perm = 5

		self.int_bnds = [[0 for x in xrange(2)] for x in xrange(5)]
		self.int_bnds[0][0] = -30
		self.int_bnds[0][1] = -10
		self.int_bnds[1][0] = -20
		self.int_bnds[1][1] = 0
		self.int_bnds[2][0] = -10
		self.int_bnds[2][1] = 10
		self.int_bnds[3][0] = 0
		self.int_bnds[3][1] = 20
		self.int_bnds[4][0] = 10
		self.int_bnds[4][1] = 30
		
		self.dbl_bnds = [[0 for x in xrange(2)] for x in xrange(5)]
		self.dbl_bnds[0][0] = -30.1
		self.dbl_bnds[0][1] = -10.2
		self.dbl_bnds[1][0] = -20.3
		self.dbl_bnds[1][1] = 0.0
		self.dbl_bnds[2][0] = -10.4
		self.dbl_bnds[2][1] = 10.5
		self.dbl_bnds[3][0] = 0.0
		self.dbl_bnds[3][1] = 20.6
		self.dbl_bnds[4][0] = 10.7
		self.dbl_bnds[4][1] = 30.8

		self.fit_sum = 0
		self.fit_avg = 0
		self.fit_max = 0
		self.fit_std_dev = 0

		self.max_chrome = GA_Chrome()
		self.chrome = [ GA_Chrome() for i in range(self.pop_size)]

	def calc_stats(self):
		max_index = -1
		self.fit_sum = 0.0

		for i in range(self.pop_size):
			self.fit_sum += self.chrome[i].fitness
			if self.chrome[i].fitness > self.fit_max:
				self.fit_max = self.chrome[i].fitness
				max_index = i

		self.fit_avg = self.fit_sum / self.pop_size

		self.fit_std_dev = 0.0
		for i in range(self.pop_size):
			temp = self.chrome[i].fitness - self.fit_avg
			temp *= temp
			self.fit_std_dev += temp

		self.fit_std_dev /= self.pop_size
		self.fit_std_dev = math.sqrt(self.fit_std_dev)

		if -1 != max_index:
			self.copy_chrome(self.chrome[max_index], self.max_chrome)

	def evolve(self):
		for i in range(self.pop_size): # create (init) population
			self.pop_chrome(self.chrome[i])
			self.chrome[i].fitness = self.evaluate(self.chrome[i])
			self.print_chrome(self.chrome[i])

		self.calc_stats()
		self.print_stats()

		for i in range(self.max_iter):
			self.parent_selection()
			# self.print_chrome(self.chrome[i])

	def fitness_proportional_selection(self):
		print 'fitness'
		# Goldberg's sigma scaling (or use windowing)
		for i in range(self.pop_size):
			temp = self.chrome[i].fitness - self.fit_avg + (self.s*self.fit_std_dev)
			self.chrome[i].fitness = temp if temp > 0.0 else 0.0
			sys.stdout.write(str(i))
			sys.stdout.write(' ')
			sys.stdout.write(str(self.chrome[i].fitness))
			sys.stdout.write(' ')
			sys.stdout.write(str(temp))
			print ' '
		print ' '

		for i in range(1, self.pop_size):
			self.chrome[i].fitness += self.chrome[i-1].fitness
			sys.stdout.write(str(self.chrome[i].fitness))
			print ' '
		print ' '

		for i in range(self.pop_size):
			self.chrome[i].fitness /= self.chrome[self.pop_size-1].fitness
			sys.stdout.write(str(self.chrome[i].fitness))
			print ' '
		print ' '

	def parent_selection(self):
		mating_pool = [ GA_Chrome() for i in range(self.pop_size)]

		# keep fittest chrome and copy it to first child chrome
		self.copy_chrome(self.max_chrome, mating_pool[0])

		# choose one of these
		if self.parent_choice == "f":
			self.fitness_proportional_selection()
		elif self.parent_choice == "r":
			pass
			# self.ranking_selection()
		else: # self.parent_choice == "t"
			pass
			# self.tournament_selection()

		# then one of these (if not "tournament")
		if self.parent_choice != "t":
			pass
			# self.selection_probability(mating_pool)

		# merge children back into population
		for i in range(1,self.pop_size):
			self.copy_chrome(mating_pool[i], self.chrome[i])

	def print_stats(self):
		print 'sum = ' + str(self.fit_sum) + \
			' avg = ' + str(self.fit_avg) + \
			' max = ' + str(self.fit_max) + \
			' std_dev = ' + str(self.fit_std_dev)
		self.print_chrome(self.max_chrome)

class GA_Chrome:
	def __init__(self):
		self.bool_chrome = []
		self.int_chrome = []
		self.dbl_chrome = []
		self.perm_chrome = []

		self.fitness = 0

	def copy_chrome(self, from_chrome, to_chrome):
		to_chrome.bool_chrome = list(from_chrome.bool_chrome)
		to_chrome.int_chrome = list(from_chrome.int_chrome)
		to_chrome.dbl_chrome = list(from_chrome.dbl_chrome)
		to_chrome.perm_chrome = list(from_chrome.perm_chrome)

		to_chrome.fitness = from_chrome.fitness

	def evaluate(self, chrome):
		fitness = 0.0

		bool_sum = 0
		for i in range(self.num_bool):
			if True == chrome.bool_chrome[i]:
				bool_sum += 1

		int_sum = 0
		for i in range(self.num_int):
			int_sum += chrome.int_chrome[i]

		dbl_sum = 0.0
		for i in range(self.num_dbl):
			dbl_sum += chrome.dbl_chrome[i]

		perm_sum = 0
		for i in range(self.num_perm):
			perm_sum += abs(chrome.perm_chrome[i] - (i+1))

		fitness = bool_sum + int_sum + dbl_sum + perm_sum

		return fitness

	def pop_chrome(self, chrome):
		for i in range(self.num_bool):
			chrome.bool_chrome.append(random.choice([True, False]))

		for i in range(self.num_int):
			chrome.int_chrome.append(random.randint(self.int_bnds[i][0], self.int_bnds[i][1]))

		for i in range(self.num_dbl):
			chrome.dbl_chrome.append(random.uniform(self.dbl_bnds[i][0], self.dbl_bnds[i][1]))

		for i in range(self.num_perm):
			chrome.perm_chrome.append(i)
		self.sattolo_cycle(chrome.perm_chrome)

	def print_chrome(self, chrome):
		for i in range(self.num_bool):
			sys.stdout.write(str(chrome.bool_chrome[i]))
			sys.stdout.write(' ')

		for i in range(self.num_int):
			sys.stdout.write(str(chrome.int_chrome[i]))
			sys.stdout.write(' ')

		for i in range(self.num_dbl):
			sys.stdout.write(str(chrome.dbl_chrome[i]))
			sys.stdout.write(' ')

		for i in range(self.num_perm):
			sys.stdout.write(str(chrome.perm_chrome[i]))
			sys.stdout.write(' ')

		sys.stdout.write(str(chrome.fitness))

		print ' '

	def sattolo_cycle(self, items): # knuth shuffle
		i = len(items)
		while i > 1:
		    i = i - 1
		    j = random.randrange(i)  # 0 <= j <= i-1
		    items[j], items[i] = items[i], items[j]

def main(args):
	ga = GA()
	ga.evolve()

if __name__ == '__main__':
	main(sys.argv)
