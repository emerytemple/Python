#!/usr/bin/python3.1

# ga
pop_size = 10
max_iter = 1

num_bool = 10
num_int = 5
num_dbl = 5
num_perm = 2

int_low_bnd = [-30, -20, -10, 0, 10]
int_high_bnd = [-10, 0, 10, 20, 30]

dbl_low_bnd = [-30.1, -20.3, -10.4, 0.0, 10.7]
dbl_high_bnd = [-10.2, 0.0, 10.5, 20.6, 30.8]

perm_len = [5,10]

# parent
parent_choice = 'ranking' # fps, ranking, tournament

# 2.0 for fps,
# 1.5 (1.0 < s <= 2.0) for ranking,
# num chosen (1 <= s <= pop_size-1) for tournament
s = 1.5

selection_prob = 'sus' # roulette, sus

# meiosis
pc = 1.0 # crossover rate

bool_cross = 'npoint' # npoint, uniform
bool_npc = 1
bool_ucp = 0.5

int_cross = 'npoint' # npoint, uniform
int_npc = 1
int_ucp = 0.5

dbl_cross = 'discrete' # discrete, simple, single, whole
alpha = 0.5

perm_cross = 'pmx' # pmx, edge, order, cycle

# chrome
bool_mutate = 'bitwise' # bitwise
bool_pm = 0.05

int_mutate = 'creep' # random, creep
int_pm = 0.5
int_sigma = 1.0

dbl_mutate = 'nonuniform' # uniform, nonuniform
dbl_pm = 0.5
dbl_sigma = 1.0

perm_mutate = 'swap' # swap, insert, scramble, inversion
perm_pm = 0.5
