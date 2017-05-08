#!/usr/bin/python3.1

import logging
import random

from math import ceil


'''	bool
		bitwise_mutation
	int
		creep_mutation
		random_resetting
	double
		uniform mutation
		nonuniform mutation
	perm
		insert mutation
		inversion mutation
		scramble mutation
		swap mutation
'''

class Chromosome:
	def __init__(self):
		logging.getLogger()

		self.bool_chrome = []
		self.int_chrome = []
		self.dbl_chrome = []
		self.perm_chrome = [[]]

		self.fitness = None

	def __str__(self):
		return str(self.bool_chrome) + "\n" + \
				str(self.int_chrome) + "\n" + \
				str(self.dbl_chrome) + "\n" + \
				str(self.perm_chrome) + "\n"

	def bitwise_mutation(self, pm):
		''' Consider each gene separately and allow each bit to flip
			with a small probability of pm.
		'''

		logging.debug('Bitwise mutation...')
		logging.debug('pm = %s', pm)

		for i, chrome in enumerate(self.bool_chrome):
			rand = random.uniform(0.0, 1.0)

			if rand < pm:
				self.bool_chrome[i] = False if chrome else True

			logging.debug('%s %s %s', chrome, rand, self.bool_chrome[i])

		logging.debug(' ')

	def creep_mutation(self, pm, ilow, ihigh, sigma):
		''' Add a small (positive or negative) value to each gene
			with probability p.
		'''

		logging.debug('Creep mutation...')
		logging.debug('pm = %s', pm)

		for i, chrome in enumerate(self.int_chrome):
			rand = random.uniform(0.0, 1.0)

			rand2 = None
			if rand < pm:
				rand2 = int(ceil(random.gauss(0.0, sigma)))
				self.int_chrome[i] += rand2

				if self.int_chrome[i] > ihigh[i]:
					self.int_chrome[i] = ihigh[i]
				elif self.int_chrome[i] < ilow[i]:
					self.int_chrome[i] = ilow[i]

			logging.debug('%s %s + %s = %s', \
				rand, chrome, rand2, self.int_chrome[i])

		logging.debug(' ')

	def evaluate(self):
		''' Calculate the fitness of a chromosome.
		'''

		bool_sum = sum([1 if i == True else 0 for i in self.bool_chrome])
		int_sum = sum(self.int_chrome)
		dbl_sum = sum(self.dbl_chrome)

		total_sum = 0
		for i in range(len(self.perm_chrome)):
			perm_sum = 0
			for j in range(len(self.perm_chrome[i])):
				perm_sum += len(self.perm_chrome[i]) - 1 \
						- abs(self.perm_chrome[i][j] - j)
			total_sum += perm_sum

		self.fitness = bool_sum + int_sum + dbl_sum + total_sum

	def insert_mutation(self, pm, perm_len):
		''' Pick two alleles at random and move one
			so that it is next to the other,
			shuffling along the others to make room.
		'''

		logging.debug('Insert mutation...')
		logging.debug('pm = %s', pm)

		for i, chrome in enumerate(self.perm_chrome):
			rand = random.uniform(0.0, 1.0)

			if rand < pm and len(chrome) > 1:
				num1, num2 = random.sample(range(perm_len[i]), 2)
				if num1 > num2:
					num1, num2 = num2, num1

				logging.debug('%s %s %s', rand, num1, num2)
				logging.debug('	%s', chrome)

				self.perm_chrome[i].insert(num1, self.perm_chrome[i].pop(num2))

				logging.debug('	%s', self.perm_chrome[i])
			else:
				logging.debug(rand)

		logging.debug(' ')

	def inversion_mutation(self, pm, perm_len):
		''' Randomly select two positions in the string
			and reverse the order in which the values appear
			between those positions.
		'''

		logging.debug('Inversion mutation...')
		logging.debug('pm = %s', pm)

		for i, chrome in enumerate(self.perm_chrome):
			rand = random.uniform(0.0, 1.0)

			if rand < pm and len(chrome) > 1:
				num1, num2 = random.sample(range(perm_len[i]), 2)
				if num1 > num2:
					num1, num2 = num2, num1

				aaa = [self.perm_chrome[i][j] for j in range(num1, num2+1)]

				logging.debug('	%s', self.perm_chrome[i])
				logging.debug('		%s', aaa)

				logging.debug('%s %s %s', rand, num1, num2)

				aaa.reverse()

				logging.debug('		%s', aaa)

				self.perm_chrome[i][num1:num2+1] = aaa

				logging.debug('	%s', self.perm_chrome[i])
			else:
				logging.debug(rand)

		logging.debug(' ')

	def nonuniform_mutation(self, pm, dlow, dhigh, sigma):
		''' Add to the current gene value an amount
			drawn randomly from a Gaussian distribution with mean zero
			and user-specified standard deviation.

			Then curtail the resulting value to the range [L, U]
			if necessary.
		'''

		logging.debug('Nonuniform mutation...')
		logging.debug('pm = %s', pm)

		for i, chrome in enumerate(self.dbl_chrome):
			rand = random.uniform(0.0, 1.0)

			rand2 = None
			if rand < pm:
				rand2 = random.gauss(0.0, sigma)
				self.dbl_chrome[i] += rand2

				if self.dbl_chrome[i] > dhigh[i]:
					self.dbl_chrome[i] = dhigh[i]
				elif self.dbl_chrome[i] < dlow[i]:
					self.dbl_chrome[i] = dlow[i]

			logging.debug('%s %s + %s = %s', \
				rand, chrome, rand2, self.dbl_chrome[i])

		logging.debug(' ')

	def pop_chrome(self, nbool, nint, ndbl, nperm, il, ih, dl, dh, pl):
		''' Initialize chromosome with random values.
		'''

		self.bool_chrome = [random.choice([True, False]) for i in range(nbool)]
		self.int_chrome = [random.randint(il[i], ih[i]) for i in range(nint)]
		self.dbl_chrome = [random.uniform(dl[i], dh[i]) for i in range(ndbl)]
		self.perm_chrome = [[j+1 for j in range(pl[i])] for i in range(nperm)]
		for i in range(nperm):
			random.shuffle(self.perm_chrome[i])

	def random_resetting(self, pm, ilow, ihigh):
		''' A new value is chosen at random from the set of permissible values
			in each position with probability pm.
		'''

		logging.debug('Random resetting...')
		logging.debug('pm = %s', pm)

		for i, chrome in enumerate(self.int_chrome):
			rand = random.uniform(0.0, 1.0)

			if rand < pm:
				self.int_chrome[i] = random.randint(ilow[i], ihigh[i])

			logging.debug('%s %s %s', chrome, rand, self.int_chrome[i])

		logging.debug(' ')

	def scramble_mutation(self, pm, perm_len):
		''' The entire string, or some randomly chosen subset
			of values within it, have their positions scrambled.
		'''

		logging.debug('Scramble mutation...')
		logging.debug('pm = %s', pm)

		for i, chrome in enumerate(self.perm_chrome):
			rand = random.uniform(0.0, 1.0)

			if rand < pm and len(chrome) > 1:
				num1, num2 = random.sample(range(perm_len[i]), 2)
				if num1 > num2:
					num1, num2 = num2, num1

				aaa = [self.perm_chrome[i][j] for j in range(num1, num2+1)]

				logging.debug('	%s', self.perm_chrome[i])
				logging.debug('		%s', aaa)

				random.shuffle(aaa)

				logging.debug('		%s', aaa)

				self.perm_chrome[i][num1:num2+1] = aaa

				logging.debug('	%s', self.perm_chrome[i])
			else:
				logging.debug(rand)

	def swap_mutation(self, pm, perm_len):
		''' Pick two random positions (genes) in the string
			and swap their allele values.
		'''

		logging.debug('Swap mutation...')
		logging.debug('pm = %s', pm)

		for i, chrome in enumerate(self.perm_chrome):
			rand = random.uniform(0.0, 1.0)

			if rand < pm and len(chrome) > 1:
				num1, num2 = random.sample(range(perm_len[i]), 2)
				if num1 > num2:
					num1, num2 = num2, num1

				logging.debug('%s %s %s', rand, num1, num2)
				logging.debug('	%s', chrome)

				self.perm_chrome[i][num1], self.perm_chrome[i][num2] = \
				self.perm_chrome[i][num2], self.perm_chrome[i][num1]

				logging.debug('	%s', self.perm_chrome[i])
			else:
				logging.debug(rand)

		logging.debug(' ')

	def uniform_mutation(self, pm, dlow, dhigh):
		''' Values are drawn randomly from[L,U].
		'''

		logging.debug('Uniform mutation...')
		logging.debug('pm = %s', pm)

		for i, chrome in enumerate(self.dbl_chrome):
			rand = random.uniform(0.0, 1.0)

			if rand < pm:
				self.dbl_chrome[i] = random.uniform(dlow[i], dhigh[i])

			logging.debug('%s %s %s', chrome, rand, self.dbl_chrome[i])

		logging.debug(' ')
