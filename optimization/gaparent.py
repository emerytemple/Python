#!/usr/bin/python3.1

import copy
import logging
import operator
import random

def fitness_proportional_selection(pop_size, chrome, fit_avg, fit_std_dev, s):
	''' Scale fitness values (with windowing) to [0,1].
		Allocate selection probabilities to individuals
		according to their actual fitness values.
	'''

	logging.debug('Fitness proportional selection...')

	log = [['' for x in range(1)] for y in range(pop_size)]

	# Goldberg's sigma scaling
	scale = fit_avg - (s*fit_std_dev)
	logging.debug('scale = %s', scale)

	arr = [max(c.fitness - scale, 0.0) for c in chrome]

	logstr = '0 ' + str(chrome[0].fitness) + ' ' \
				  + str(arr[0]) + ' ' + str(arr[0]) + ' '
	logging.debug(logstr)
	for i in range(1, pop_size):
		logstr = str(i) + ' ' + str(chrome[i].fitness) + ' ' + str(arr[i]) + ' '
		arr[i] += arr[i-1]
		logstr += str(arr[i])
		logging.debug(logstr)
	logging.debug(' ')

	max_val = arr[pop_size-1]
	arr3 = [i/max_val for i in arr]
	for i in range(pop_size):
		chrome[i].fitness = arr3[i]
		logging.debug(str(i) + ' ' + str(arr3[i]))

	logging.debug(' ')

def ranking_selection(pop_size, chrome, s):
	''' Allocate selection probabilities to individuals
		according to their rank.
	'''

	logging.debug('Ranking selection...')

	log = [['' for x in range(1)] for y in range(pop_size)]
	temp = [[0 for x in range(2)] for y in range(pop_size)]

	for i in range(pop_size):
		temp[i][0] = i
		temp[i][1] = chrome[i].fitness

		log[i].append(str(temp[i][0]))
		log[i].append(' ')
		log[i].append(str(temp[i][1]))
		log[i].append(' ')

	# sort based on fitness
	temp = sorted(temp, key=operator.itemgetter(1))

	for i in range(pop_size):
		log[i].append(str(temp[i][0]))
		log[i].append(' ')
		log[i].append(str(temp[i][1]))
		log[i].append(' ')

	# calculate selection probability
	part1 = (2.0-s)/pop_size
	part2 = (2.0*(s-1.0))
	part3 = pop_size*(pop_size-1.0)
	for i in range(pop_size):
		temp[i][1] = part1 + ((i*part2)/part3)

		log[i].append(str(temp[i][0]))
		log[i].append(' ')
		log[i].append(str(temp[i][1]))
		log[i].append(' ')

	# sort based on chrome number
	temp = sorted(temp, key=operator.itemgetter(0))

	for i in range(pop_size):
		log[i].append(str(temp[i][0]))
		log[i].append(' ')
		log[i].append(str(temp[i][1]))
		log[i].append(' ')

	# normalize fitness array
	log[0].append('0 ')
	log[0].append(str(temp[0][1]))
	for i in range(1,pop_size):
		temp[i][1] += temp[i-1][1]

		log[i].append(str(temp[i][0]))
		log[i].append(' ')
		log[i].append(str(temp[i][1]))
		log[i].append(' ')

	for i in range(pop_size):
		chrome[i].fitness = temp[i][1]

	for i in range(pop_size):
		logcat = ''.join(log[i])
		logging.debug(logcat)

	logging.debug(' ')

def roulette_wheel_sampling(pop_size, chrome, mating_pool):
	''' Spin a one-armed roulette wheel, 
		where the sizes of the holes reflect the selection probabilities.
	'''

	logging.debug('Roulette wheel algorithm...')

	log = [['' for x in range(1)] for y in range(pop_size)]

	for i in range(pop_size):
		log[i].append(str(i))
		log[i].append(' ')

		rndm = random.uniform(0.0, 1.0)

		log[i].append(str(rndm))
		log[i].append(' ')

		temp = 0
		for j in range(pop_size):
			if rndm > chrome[j].fitness:
				temp = j+1

		log[i].append(str(temp))

		mating_pool[i] = copy.deepcopy(chrome[temp])

	for i in range(pop_size):
		logcat = ''.join(log[i])
		logging.debug(logcat)

	logging.debug(' ')

def stochastic_universal_sampling(pop_size, chrome, mating_pool):
	''' Make one spin of a wheel with equally spaced arms.
		Number of arms is pop_size.
	'''

	log = [['' for x in range(1)] for y in range(pop_size)]
	logging.debug('Stochastic universal sampling...')

	for i in range(pop_size):
		log[i].append(str(i))
		log[i].append(' ')

		if 0 == i:
			rndm = random.uniform(0.0, 1.0)/(pop_size)
		else:
			rndm += 1.0/(pop_size)

		log[i].append(str(rndm))
		log[i].append(' ')

		temp = 0
		for j in range(pop_size):
			if rndm > chrome[j].fitness:
				temp = j+1

		log[i].append(str(temp))

		mating_pool[i] = copy.deepcopy(chrome[temp])

	for i in range(pop_size):
		logcat = ''.join(log[i])
		logging.debug(logcat)

	logging.debug(' ')

def tournament_selection(pop_size, chrome, mating_pool, tsk):
	''' Subjectively select among the individuals
		by comparing their fitness values
		without requiring any global knowledge of the population.

		TODO: add replacement flag
	'''

	logging.debug('Tournament selection...')

	log = [['' for x in range(1)] for y in range(pop_size)]

	for i in range(pop_size):
		log[i].append(str(i))
		log[i].append(' ')

		chosen_index = [random.randint(0, pop_size-1) for j in range(tsk)]
		log[i].append(str(chosen_index))
		log[i].append(' ')

		chosen_fitness = [chrome[j].fitness for j in chosen_index]
		log[i].append(str(chosen_fitness))
		log[i].append(' ')

		winner = chosen_fitness.index(max(chosen_fitness))
		log[i].append(str(winner))
		log[i].append(' ')
		log[i].append(str(chosen_index[winner]))

		mating_pool[i] = copy.deepcopy(chrome[chosen_index[winner]])

	for i in range(pop_size):
		logcat = ''.join(log[i])
		logging.debug(logcat)

	logging.debug(' ')
