#!/usr/bin/python3.1

import copy
import logging
import math
import random
# import stat

import gachrome as gac
import gameiosis as gam
import gaparent as gap
import gaparam as param

# import matplotlib.pyplot as plt

class GA: # population
	def __init__(self):
		logging.getLogger()

		self.fit_sum = 0
		self.fit_avg = 0
		self.fit_max = 0
		self.fit_std_dev = 0

		self.max_chrome = gac.Chromosome()
		self.chrome = [gac.Chromosome() for i in range(param.pop_size)]

	def calc_stats(self):
		''' Calculate sum, avg, max, std_dev, of entire population
			and index of most fit chromesome.
		'''

		logging.debug('Calculating statistics of entire population...')

		arr = [c.fitness for c in self.chrome]
		logging.debug(arr)

		self.fit_sum = sum(arr)
		self.fit_avg = self.fit_sum/param.pop_size
		self.fit_max = max(arr)

		''' use built-in statistics module for stdev()
			after upgrading from ubuntu 10.04
			so you can upgrade to python 3.4
		'''
		mysum = 0.0
		for fitness in arr:
			mysum += (fitness-self.fit_avg)**2
		self.fit_std_dev = math.sqrt(mysum/len(arr))

		logging.info('sum = %s, avg = %s, max = %s, std_dev = %s', \
			self.fit_sum, self.fit_avg, self.fit_max, self.fit_std_dev)
		logging.debug(' ')

		return self.fit_sum, self.fit_avg, self.fit_max, self.fit_std_dev

	def elitism(self):
		''' Prevent the loss of the current fitest member
			of the population. (keep a trace of the current fittest member)
		'''

		logging.debug('Elitism...')

		arr = [c.fitness for c in self.chrome]
		logging.debug(arr)

		max_index = arr.index(max(arr))
		logging.debug('max_index = %s', max_index)

		if self.max_chrome.fitness == None:
			logging.debug('Finding initial max_chrome...')
			self.max_chrome = copy.deepcopy(self.chrome[max_index])
			logging.debug(self.max_chrome.__str__())
		elif self.chrome[max_index].fitness >= self.max_chrome.fitness:
			logging.debug('Saving new max_chrome...')
			self.max_chrome = copy.deepcopy(self.chrome[max_index])
			logging.debug(self.max_chrome.__str__())
		else:
			logging.debug('No new max, discarding offspring to make room...')
			self.chrome[0] = copy.deepcopy(self.max_chrome)

		logging.debug(' ')

	def evolve(self):
		''' Create new generations of chromosomes. (main loop for GA)
		'''

		fsum = [0 for i in range(param.max_iter+1)]
		favg = [0 for i in range(param.max_iter+1)]
		fmax = [0 for i in range(param.max_iter+1)]
		fsdv = [0 for i in range(param.max_iter+1)]

		xxx = [i for i in range(param.max_iter+1)]

		self.initialize()
		self.elitism()
		fsum[0], favg[0], fmax[0], fsdv[0] = self.calc_stats()

		for i in range(param.max_iter):
			logging.debug('Iteration Number %s...', i)
			logging.debug(' ')

			self.parent_selection()
			self.recombination()
			self.mutation()

			for j in range(param.pop_size):
				self.chrome[j].evaluate()

			self.elitism()
			fsum[i+1], favg[i+1], fmax[i+1], fsdv[i+1] = self.calc_stats()

		print(self.max_chrome)

		'''
		plt.title('Genetic Algorithm Progress')
		plt.xlabel('Iteration')
		plt.ylabel('Performance')
		p1 = plt.plot(xxx, favg)
		p2 = plt.plot(xxx, fmax)
		plt.legend([p1, p2], ["avg", "fmax"], loc='lower right')
		plt.show()
		'''

	def initialize(self):
		''' Populate the chromes and find their fitness.
		'''

		logging.debug('Creating initial population...')

		for i in range(param.pop_size):
			self.chrome[i].pop_chrome(param.num_bool, \
									  param.num_int, \
									  param.num_dbl, \
									  param.num_perm, \
									  param.int_low_bnd, \
									  param.int_high_bnd, \
									  param.dbl_low_bnd, \
									  param.dbl_high_bnd, \
									  param.perm_len)
			self.chrome[i].evaluate()
			logging.debug(self.chrome[i].__str__())

		logging.debug(' ')

	def mutation(self):
		''' Use only one parent and create one child
			by applying some kind of randomised change
			to the representation.
		'''

		logging.debug('Mutation...')

		for i in range(param.pop_size):
			# bool
			if 0 != param.num_bool:
				logging.debug('Bool mutation...')
				self.chrome[i].bitwise_mutation(param.bool_pm)

			# int
			if 0 != param.num_int:
				if param.int_mutate == 'random':
					self.chrome[i].random_resetting(param.int_pm, \
													param.int_low_bnd, \
													param.int_high_bnd)
				elif param.int_mutate == 'creep':
					self.chrome[i].creep_mutation(param.int_pm, \
												  param.int_low_bnd, \
												  param.int_high_bnd, \
												  param.int_sigma)
				else:
					pass

			# double
			if 0 != param.num_dbl:
				if param.dbl_mutate == 'uniform':
					self.chrome[i].uniform_mutation(param.dbl_pm, \
													param.dbl_low_bnd, \
													param.dbl_high_bnd)
				elif param.dbl_mutate == 'nonuniform':
					self.chrome[i].nonuniform_mutation(param.dbl_pm, \
													   param.dbl_low_bnd, \
													   param.dbl_high_bnd, \
													   param.dbl_sigma)
				else:
					pass

			# perm
			if 0 != param.num_perm:
				if param.perm_mutate == 'swap':
					self.chrome[i].swap_mutation(param.perm_pm, \
												 param.perm_len)
				elif param.perm_mutate == 'insert':
					self.chrome[i].insert_mutation(param.perm_pm, \
												   param.perm_len)
				elif param.perm_mutate == 'scramble':
					self.chrome[i].scramble_mutation(param.perm_pm, \
													 param.perm_len)
				elif param.perm_mutate == 'inversion':
					self.chrome[i].inversion_mutation(param.perm_pm, \
													  param.perm_len)
				else:
					pass

	def parent_selection(self):
		''' Select parents to create offspring from existing population.
		'''

		logging.debug('Selecting parents for breeding...')

		mating_pool = [ gac.Chromosome() for i in range(param.pop_size)]

		# populate rest of chromes
		if param.parent_choice == 'fps':
			gap.fitness_proportional_selection(param.pop_size, \
											   self.chrome, \
											   self.fit_avg, \
											   self.fit_std_dev, \
											   param.s)
		elif param.parent_choice == 'ranking':
			gap.ranking_selection(param.pop_size, \
								  self.chrome, \
								  param.s)
		elif param.parent_choice == 'tournament':
			gap.tournament_selection(param.pop_size, \
									 self.chrome, \
									 mating_pool, \
									 param.s)
		else:
			pass

		if param.parent_choice != 'tournament':
			if param.selection_prob == 'roulette':
				gap.roulette_wheel_sampling(param.pop_size, \
											self.chrome, \
											mating_pool)
			elif param.selection_prob == 'sus':
				gap.stochastic_universal_sampling(param.pop_size, \
												  self.chrome, \
												  mating_pool)
			else:
				pass

		logging.debug('Mating pool...')
		for i in range(param.pop_size):
			logging.debug(str(mating_pool[i]))
		logging.debug(' ')

		# replace entire generation by its offspring
		for i in range(param.pop_size):
			self.chrome[i] = copy.deepcopy(mating_pool[i])

	def recombination(self):
		''' Create new individuals from information contained within two
			(or more) parents.
		'''

		logging.debug('Recombination...')

		if 0 != param.num_int:
			max_int_bnds = gam.calc_max_bnds(param.num_int, \
											 param.int_low_bnd, \
											 param.int_high_bnd)

		for i in range(0,param.pop_size-1,2):
			rndm = random.uniform(0.0, 1.0)
			logging.debug('%s %s %s', i, i+1, rndm)
			logging.debug(' ')

			if rndm <= param.pc:
				# bool
				if 0 != param.num_bool:
					logging.debug('Bool crossover...')

					if param.bool_cross == 'npoint':
						gam.npoint_crossover(self.chrome[i].bool_chrome, \
											 self.chrome[i+1].bool_chrome, \
											 param.bool_npc)
					elif param.bool_cross == 'uniform':
						gam.uniform_crossover(self.chrome[i].bool_chrome, \
											 self.chrome[i+1].bool_chrome, \
											 param.bool_ucp)
					else:
						pass

				# int
				if 0 != param.num_int:
					logging.debug('Int crossover...')
				for j in range(param.num_int):
					chrome1, chrome2, mag1, mag2 = gam.encode(max_int_bnds[j], \
												self.chrome[i].int_chrome[j], \
												self.chrome[i+1].int_chrome[j])

					if param.int_cross == 'npoint':
						gam.npoint_crossover(chrome1, chrome2, param.int_npc)
					elif param.int_cross == 'uniform':
						gam.uniform_crossover(chrome1, chrome2, param.int_ucp)
					else:
						pass

					ccc1 = gam.decode(chrome1, mag1)
					ccc2 = gam.decode(chrome2, mag2)

					if param.int_low_bnd[j] <= ccc1 \
						and ccc1 <= param.int_high_bnd[j]:
						self.chrome[i].int_chrome[j] = ccc1

					if param.int_low_bnd[j] <= ccc2 \
						and ccc2 <= param.int_high_bnd[j]:
						self.chrome[i+1].int_chrome[j] = ccc2

				# double
				if 0 != param.num_dbl:
					logging.debug('Double crossover...')

					chrome1 = self.chrome[i].dbl_chrome
					chrome2 = self.chrome[i+1].dbl_chrome

					if param.dbl_cross == 'discrete':
						gam.discrete_recombination(chrome1, chrome2)
					elif param.dbl_cross == 'simple':
						gam.simple_recombination(chrome1, chrome2, param.alpha)
					elif param.dbl_cross == 'single':
						gam.single_recombination(chrome1, chrome2, param.alpha)
					elif param.dbl_cross == 'whole':
						gam.whole_recombination(chrome1, chrome2, param.alpha)
					else:
						pass

					self.chrome[i].dbl_chrome = chrome1
					self.chrome[i+1].dbl_chrome = chrome2

				# perm
				if 0 != param.num_perm:
					logging.debug('Perm crossover...')
				for j in range(param.num_perm):
					chrome1 = self.chrome[i].perm_chrome[j]
					chrome2 = self.chrome[i+1].perm_chrome[j]

					if param.perm_cross == 'pmx':
						gam.partially_mapped_crossover(chrome1, chrome2)
					elif param.perm_cross == 'edge':
						gam.edge_crossover(chrome1, chrome2)
					elif param.perm_cross == 'order':
						gam.order_crossover(chrome1, chrome2)
					elif param.perm_cross == 'cycle':
						gam.cycle_crossover(chrome1, chrome2)
					else:
						pass
