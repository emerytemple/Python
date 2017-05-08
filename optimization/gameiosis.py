#!/usr/bin/python3.1

import logging
import random

import gachrome as gac

def calc_max_bnds(num_int, int_low_bnd, int_high_bnd):
	''' Convert int_bnds to binary, find the max of the lengths
		of the binary numbers. Do this for each gene.
	'''

	max_int_bnds = []
	for i in range(num_int):
		mag1 = 1 if int_low_bnd[i] > 0 else 0
		bc1 = bin(int_low_bnd[i])[2 if mag1 > 0 else 3:]
		len1 = len(bc1)

		mag2 = 1 if int_high_bnd[i] > 0 else 0
		bc2 = bin(int_high_bnd[i])[2 if mag2 > 0 else 3:]
		len2 = len(bc2)

		max_len = len1 if len1 > len2 else len2
		max_int_bnds.append(max_len)

	return max_int_bnds

def cycle_crossover(chrome1, chrome2):
	''' Preserve as much information as possible about the
		absolute position in which elements occur.
		Divide the elements into cycles. Then create offspring
		by selecting alternate cycles from each parent.
	'''

	logging.debug('Cycle crossover...')

	length = len(chrome1)

	c1 = [-1 for i in range(length)]

	logging.debug('%s', chrome1)
	logging.debug('%s', chrome2)
	logging.debug(' ')

	inst = 0
	for i in range(length):
		logging.debug('i = %s, c1[%s] = %s', i, i, c1[i])
		if c1[i] == -1:
			temp = i
			c1[i] = inst
			logging.debug('	c1[%s] = %s', i, c1[i])
			for j in range(length):
				ind = chrome1.index(chrome2[temp])
				logging.debug('		j = %s	start = %s	\
									from = %s	to = %s	ind = %s', \
									j, chrome1[temp], chrome2[temp], \
									chrome1[chrome1.index(chrome2[temp])], \
									chrome1.index(chrome2[temp]))
				c1[ind] = inst
				if ind == i:
					inst += 1
					break
				else:
					temp = ind

	logging.debug('%s	%s', chrome1, c1)
	logging.debug('%s', chrome2)
	logging.debug(' ')

	for i in range(length):
		if c1[i] % 2 != 0:
			chrome1[i], chrome2[i] = chrome2[i], chrome1[i]

	logging.debug('%s', chrome1)
	logging.debug('%s', chrome2)
	logging.debug(' ')

def decode(chrome, mag):
	''' Convert binary string to integer value,
		changing sign if necessary.
	'''

	logging.debug('Decode...')

	cf = int(''.join(['1' if x == True else '0' for x in chrome]), base=2)
	ccc = cf if mag == 1 else -cf

	return ccc

def discrete_recombination(chrome1, chrome2):
	''' An allele is one floating-point value instead of one bit.
	'''

	logging.debug('Discrete recombination...')

	length = len(chrome1)
	rand = random.randint(1, length-1)
	logging.debug('rand = %s', rand)

	logging.debug(chrome1)
	logging.debug(chrome2)
	logging.debug(' ')

	for j in range(rand,length):
		chrome1[j], chrome2[j] = chrome2[j], chrome1[j]

	logging.debug(chrome1)
	logging.debug(chrome2)

	logging.debug(' ')

def edge_crossover(chrome1, chrome2):
	''' Create offspring using only edges
		that are present in one or more parent.
	'''

	logging.debug('Edge crossover...')

	# construct edge table
	length = len(chrome1)

	tab = [[] for i in chrome1]

	for i in range(length):
		loc = chrome1.index(i)
		if 0 == loc:
			ind1 = length-1
			ind2 = loc+1
		elif length-1 == loc:
			ind1 = loc-1
			ind2 = 0
		else:
			ind1 = loc-1
			ind2 = loc+1
		tab[i].append(chrome1[ind1])
		tab[i].append(chrome1[ind2])

		loc = chrome2.index(i)
		if 0 == loc:
			ind1 = length-1
			ind2 = loc+1
		elif length-1 == loc:
			ind1 = loc-1
			ind2 = 0
		else:
			ind1 = loc-1
			ind2 = loc+1

		tab[i].append(chrome2[ind1])
		tab[i].append(chrome2[ind2])

	# create offspring
	offspring = []
	for i in range(length):
		if 0 == i: # pick initial random element
			curr = random.randint(0, length-1)
		else:
			lenarr = []
			for j in range(len(tab[curr])):
				aaa = list(set(tab[tab[curr][j]]))
				lenarr.append(len(aaa))

			if len(tab[curr]) == 1:
				curr = tab[curr][0]
			elif len(tab[curr]) != len(list(set(tab[curr]))):
				ddd = list(set([x for x in tab[curr] \
								if tab[curr].count(x)  > 1]))
				curr = random.choice(ddd)
			else:
				slist = [x for x in tab[curr] \
						if lenarr[tab[curr].index(x)] == min(lenarr)]
				selectlist = []
				for k in range(len(lenarr)):
					if lenarr[k] == min(lenarr):
						selectlist.append(tab[curr][k])
				curr = random.choice(selectlist)

		offspring.append(curr)

		for k in range(length):
			tab[k] = [j for j in tab[k] if j != curr]

	chrome1 = offspring
	chrome2 = offspring

	logging.debug(' ')

def encode(max_len, chrome1, chrome2):
	''' Convert two integers to binary strings of the same length
		and keep track of their magnitude.
	'''

	logging.debug('Encode...')

	mag1 = 1 if chrome1 > 0 else 0
	bc1 = bin(chrome1)[2 if mag1 > 0 else 3:]
	bcf1 = bc1.zfill(max_len)
	bcs1 = [int(x) for x in bc1.zfill(max_len)]
	c1 = [True if x == 1 else False for x in bcs1]

	mag2 = 1 if chrome2 > 0 else 0
	bc2 = bin(chrome2)[2 if mag2 > 0 else 3:]
	bcf2 = bc2.zfill(max_len)				
	bcs2 = [int(x) for x in bc2.zfill(max_len)]
	c2 = [True if x == 1 else False for x in bcs2]

	return c1, c2, mag1, mag2

def npoint_crossover(chrome1, chrome2, np):
	''' Divide the parents into a number of sections of continuous genes
		and reassemble them to produce offspring.

		TODO: can conceivably take a long time to converge here 
		if there are a large number of random values possible
		and bool_npc is close to num_bool. Optimize algorithm later.
		(try using knuth shuffle example
		for removal of random number from list)
	'''

	logging.debug('%s-point crossover', np)

	visited = [False for i in range(len(chrome1))]

	logging.debug([1 if chrome else 0 for chrome in chrome1])
	logging.debug([1 if chrome else 0 for chrome in chrome2])
	logging.debug(' ')

	length = len(chrome1)
	while sum(visited) < np:
		rand = random.randint(1, length-1)

		if False == visited[rand]:
			visited[rand] = True

			for j in range(rand,length):
				chrome1[j], chrome2[j] = chrome2[j], chrome1[j]

	logging.debug([1 if visit else 0 for visit in visited])
	logging.debug(' ')
	logging.debug([1 if chrome else 0 for chrome in chrome1])
	logging.debug([1 if chrome else 0 for chrome in chrome2])		

	logging.debug(' ')

def order_crossover(chrome1, chrome2):
	'''
	'''

	logging.debug('Order crossover...')

	# step 1
	length = len(chrome1)

	c1 = [-1 for i in range(length)]
	c2 = [-1 for i in range(length)]

	rand1, rand2 = random.sample(range(1,length), 2)
	if rand1 > rand2:
		rand1, rand2 = rand2, rand1

	for i in range(rand1,rand2):
		c1[i] = chrome1[i]
		c2[i] = chrome2[i]

	logging.debug('%s	%s', rand1, rand2)
	logging.debug('%s	%s', chrome1, c1)
	logging.debug('%s	%s', chrome2, c2)
	logging.debug(' ')

	# step 2
	part1 = [chrome1[i] for i in range(rand1, rand2)]
	for i in range(length):
		if i == 0:
			ind = rand2
			temp = rand2
		elif temp == length-1:
			temp = 0
		else:
			temp += 1
		logging.debug('i = %s, temp = %s', i, temp)

		if part1.count(chrome2[temp]) == 0:
			logging.debug('	ch1[%s] = %s, ch2[%s] = %s, nind = %s', \
							ind, chrome1[ind], temp, chrome2[temp], ind+1)
			c1[ind] = chrome2[temp]
			if ind == length-1:
				ind = 0
			else:
				ind += 1

	# step 2.5
	part1 = [chrome2[i] for i in range(rand1, rand2)]
	for i in range(length):
		if i == 0:
			ind = rand2
			temp = rand2
		elif temp == length-1:
			temp = 0
		else:
			temp += 1
		logging.debug('i = %s, temp = %s', i, temp)

		if part1.count(chrome1[temp]) == 0:
			logging.debug('	ch1[%s] = %s, ch2[%s] = %s, nind = %s', \
							ind, chrome2[ind], temp, chrome1[temp], ind+1)
			c2[ind] = chrome1[temp]
			if ind == length-1:
				ind = 0
			else:
				ind += 1

	logging.debug('%s	%s', chrome1, c1)
	logging.debug('%s	%s', chrome2, c2)
	logging.debug(' ')

	chrome1, chrome2 = c1, c2

	logging.debug(chrome1)
	logging.debug(chrome2)
	logging.debug(' ')

def partially_mapped_crossover(chrome1, chrome2):
	'''
	'''

	logging.debug('Partially mapped crossover...')

	# step 1
	length = len(chrome1)

	c1 = [-1 for i in range(length)]
	c2 = [-1 for i in range(length)]

	rand1, rand2 = random.sample(range(1,length), 2)
	if rand1 > rand2:
		rand1, rand2 = rand2, rand1

	for i in range(rand1,rand2):
		c1[i] = chrome1[i]
		c2[i] = chrome2[i]

	logging.debug('%s	%s', rand1, rand2)
	logging.debug('%s	%s', chrome1, c1)
	logging.debug('%s	%s', chrome2, c2)
	logging.debug(' ')

	# step 2
	part1 = [chrome1[i] for i in range(rand1, rand2)]
	for i in range(rand1,rand2):
		logging.debug('i = %s	%s isIn %s = %s', \
					   i, chrome2[i], part1, part1.count(chrome2[i]))
		if part1.count(chrome2[i]) == 0:
			temp = i
			for j in range(rand1, rand2):
				ind = chrome2.index(chrome1[temp])
				logging.debug( \
					'j = %s	start = %s	from = %s	to = %s	ind = %s', \
					 j, chrome2[temp], chrome1[temp], \
					 chrome2[chrome2.index(chrome1[temp])], \
				 	 chrome2.index(chrome1[temp]))
				if ind < rand1 or rand2 <= ind:
					break
				else:
					temp = ind
			c1[ind] = chrome2[i]
			logging.debug('c1[%s] = %s', ind, chrome2[i])

	logging.debug(c1)
	logging.debug(' ')

	# step 2.5
	part2 = [chrome2[i] for i in range(rand1, rand2)]
	for i in range(rand1,rand2):
		logging.debug('i = %s	%s isIn %s = %s', \
					   i, chrome1[i], part2, part2.count(chrome1[i]))
		if part2.count(chrome1[i]) == 0:
			temp = i
			for j in range(rand1, rand2):
				ind = chrome1.index(chrome2[temp])
				logging.debug( \
					'j = %s	start = %s	from = %s	to = %s	ind = %s', \
					 j, chrome1[temp], chrome2[temp], \
					 chrome1[chrome1.index(chrome2[temp])], \
					 chrome1.index(chrome2[temp]))
				if ind < rand1 or rand2 <= ind:
					break
				else:
					temp = ind
			c2[ind] = chrome1[i]
			logging.debug('c1[%s] = %s', ind, chrome1[i])

	logging.debug(c2)
	logging.debug(' ')

	# step 3
	for i in range(length):
		if -1 == c1[i]:
			c1[i] = chrome2[i]
		if -1 == c2[i]:
			c2[i] = chrome1[i]

	logging.debug(c1)
	logging.debug(c2)

	logging.debug(' ')

	chrome1 = c1
	chrome2 = c2

def simple_recombination(chrome1, chrome2, alpha):
	''' Swap after a certain point.
	'''

	logging.debug('Simple recombination...')

	length = len(chrome1)
	rand = random.randint(1, length-1)
	logging.debug('rand = %s', rand)

	logging.debug(chrome1)
	logging.debug(chrome2)
	logging.debug(' ')

	for j in range(rand,length):
		chrome1[j] = ((1.0-alpha)*chrome1[j]) + (alpha*chrome2[j])
		chrome2[j] = (alpha*chrome1[j]) + ((1.0-alpha)*chrome2[j])

	logging.debug(chrome1)
	logging.debug(chrome2)

	logging.debug(' ')

def single_recombination(chrome1, chrome2, alpha):
	''' Swap at a random point.
	'''

	logging.debug('Single recombination...')

	rand = random.randint(1, len(chrome1)-1)
	logging.debug('rand = %s', rand)

	logging.debug(chrome1)
	logging.debug(chrome2)
	logging.debug(' ')

	chrome1[rand] = ((1.0-alpha)*chrome1[rand]) + (alpha*chrome2[rand])
	chrome2[rand] = (alpha*chrome1[rand]) + ((1.0-alpha)*chrome2[rand])

	logging.debug(chrome1)
	logging.debug(chrome2)

	logging.debug(' ')

def uniform_crossover(chrome1, chrome2, pm):
	''' Treat each gene independently and make a random choice 
		as to which parent it should be inherited from.
		In each position, if the value is below bool_ucp,
		the gene is inherited from the first parent;
		otherwise from the second.
	'''

	logging.debug('Uniform crossover...')
	logging.debug(pm)

	for i, _ in enumerate(chrome1):
		rand = random.uniform(0.0, 1.0)

		logstr = str(chrome1[i]) + ' ' + str(chrome2[i]) + ' '

		if rand < pm:
			chrome1[i], chrome2[i] = chrome2[i], chrome1[i]

		logstr += str(rand) + ' ' + str(chrome1[i]) + ' ' + str(chrome2[i])
		logging.debug(logstr)

	logging.debug(' ')

def whole_recombination(chrome1, chrome2, alpha):
	''' Swap all.
	'''

	logging.debug('Whole recombination...')
	logging.debug('alpha = %s', alpha)

	logging.debug(chrome1)
	logging.debug(chrome2)
	logging.debug(' ')

	for i, _ in enumerate(chrome1):
		chrome1[i] = ((1.0-alpha)*chrome1[i]) + (alpha*chrome2[i])
		chrome2[i] = (alpha*chrome1[i]) + ((1.0-alpha)*chrome2[i])

	logging.debug(chrome1)
	logging.debug(chrome2)

	logging.debug(' ')
