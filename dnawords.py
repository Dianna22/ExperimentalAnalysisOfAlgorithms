import cProfile, pstats
from io import StringIO
import numpy as np
import matplotlib.pyplot as plt
import itertools as it
import sys
import time


dna_codes = ['A', 'T', 'C', 'G']

DEFAULT_D = 2 # TODO D
# n = length of words
# D - hemming distance
# find the number of words of the solution


# TODO PAG 7

def hd(w1,w2):
	try:
		return sum([c1!=c2 for c1,c2 in zip(w1,w2)])
	except:
		import ipdb
		ipdb.set_trace(context=10)
		return 0

def check_constraint_hd(S, default_d=DEFAULT_D):
	"""
	(T/F, count)
	"""
	result = True
	count = 0
	for i in range(len(S)):
		for j in range(i+1, len(S)):
			if hd(S[i], S[j]) < default_d:
				count += 1
				result = False
	return (result, count)


def percent(word):
	"""
	check if exactly 50% of the nucleotides are C/G
	"""
	count = 0
	total = len(word)
	half = total/2 if total%2==0 else round(total/2)-1
	for ch in word:
		if ch in ['C', 'G']:
			count +=1
	if count == half:
		return True
	return False


def check_contraint_percent(S, default_d=None):
	"""
	percent assumend 50%
	"""
	result, count = True, 0
	for word in S:
		if not percent(word):
			result = False
			count += 1
	return (result, count)


def complement_ch(ch):
	return dna_codes[(dna_codes.index(ch) + 2) % 4]


def complement_w(w):
	return [complement_ch(ch) for ch in w]


def check_constraint_chd(S, default_d=DEFAULT_D):
	result = True
	count = 0
	for i in range(len(S)):
		for j in range(i+1, len(S)):
			if hd(S[i], complement_w(S[j])) < default_d:
				count += 1
				result = False
	return (result, count)


def check_constraints(S):
	return check_constraint_hd(S)[0] and check_contraint_percent(S)[0] and check_constraint_chd(S)[0]


def count_conflict_constraints(S, constraints=[check_constraint_hd, check_contraint_percent, check_constraint_chd]):
	count = 0
	for constraint in constraints:
		count += constraint(S)[1]
	return count


def all_others(c):
	return [dna_codes[i] for i in range(len(dna_codes)) if dna_codes[i] != c]


def all_substitutes_one_char(w):
	res = []
	for i in range(len(w)):
		for x in all_others(w[i]):
			copy_w = list(w)
			copy_w[i] = x
			res.append(copy_w)
	return res

# TOOL PT STABILIREA PARAMETRILOR (MAX_TRIES AND MAX_STEPS)
# LINK IN LAB
def stochastic_local_search(k,n):
	"""
	Generate seq of k words of length n
	Constraints:
		HD(w1,w2)>=d
		50% is C/G
		HD(w1,wcc(w2))>=d (Watson-Crick complement ~ A-T, C-G)
	"""
	# import ipdb
	# ipdb.set_trace(context=10)
	max_tries, max_steps = 30, 100
	probability = 0.3
	code_pool_length = len(dna_codes)
	best_S = []
	for i in range(max_tries):
		S = [[dna_codes[np.random.randint(code_pool_length)] for x in range(n)] for x in range(k)]
		best_S = list(S)
		for j in range(max_steps):
			if check_constraints(S):
				return (S, count_conflict_constraints(S))
			not_sat_constraints = count_conflict_constraints(S)
			# randomly select 2 words that violate one of the constraints
			index1, index2=np.random.randint(k), np.random.randint(k)
			w1, w2 = S[index1], S[index2]
			while check_constraints([w1,w2]):
				index1, index2=np.random.randint(k), np.random.randint(k)
				w1, w2 = S[index1], S[index2]
			# all words by subst one base
			M1, M2 = all_substitutes_one_char(w1), all_substitutes_one_char(w2)
			Munion = M1 + M2
			w = None
			if np.random.uniform() < probability:
				temp = np.random.randint(len(Munion))
				w = Munion[temp]
			else:
				max_violations_resolved = 0
				for word in Munion:
					copy_S = list(S)
					if word in M1:
						copy_S[index1] = word
					else:
						copy_S[index2] = word
					new_sat = count_conflict_constraints(copy_S)
					if not_sat_constraints - new_sat >= max_violations_resolved:
						max_violations_resolved = new_sat
						w = list(word)
			assert w is not None
			if w in M1:
				S[index1] = w
			else:
				S[index2] = w
			if count_conflict_constraints(S) < count_conflict_constraints(best_S):
				best_S = list(S)
	return (best_S, count_conflict_constraints(best_S))


def add_best_word(S, n):
	words = list(it.combinations_with_replacement(dna_codes, n))
	min_conflicts = sys.maxsize
	solution = None
	for word in words:
		temp = list(S)
		temp.append(list(word))
		conflicts = count_conflict_constraints(temp)
		if conflicts < min_conflicts:
			min_conflicts = conflicts
			solution = list(temp)
	assert solution is not None
	return solution

# STOP BOTH AFTER 2 MINS

# sGREEDY START WITH THE MAX LENGTH, AND THEN OPTIMIZE
def stochastic_greedy(k,n):
	"""
	Generate seq of k words of length n
	Constraints:
		HD(w1,w2)>=d
		50% is C/G
		HD(w1,wcc(w2))>=d (Watson-Crick complement ~ A-T, C-G)
	"""
	dna_codes_len = len(dna_codes)
	S = [[dna_codes[np.random.randint(dna_codes_len)] for i in range(n)]]
	for i in range(1, k):
		S = add_best_word(S,n)
	return [S, count_conflict_constraints(S)]


def plot():
	times = []
	for i in range(1,11):
		temps = []
		for j in range(10):
			t1 =time.time()
			# stochastic_local_search(i,5)
			stochastic_greedy(i,5)
			t2 = time.time()
			temps.append(t2-t1)
		times.append(np.mean(temps))
	plt.plot(list(range(1,11)), times)
	plt.plot(list(range(1,11)), times, 'ro')
	plt.title("Avg exe time of 10 runs of Stochastic Greedy")
	plt.show()

# plot()
def profiling(function_name):
	for i in range(1,11):
		stats_file = 'profiling_stats%s' %str(i)
		# cProfile.run('%s(%s,5)'%(function_name, str(i)), stats_file)
		cProfile.run('%s(%s,5)'%(function_name, str(i)))
		ps = pstats.Stats(stats_file)
		ps.sort_stats('cumulative').print_stats(1)

profiling('stochastic_greedy')
print("#####################")
profiling('stochastic_local_search')

def find_best_solutions():
	ls_avg_times, gr_avg_times = [], []
	ls_avg_sol, gr_avg_sol = [], []
	for i in range(1,11):
		ls_times, ls_best_sol = [], []
		gr_times, gr_best_sol = [], []
		for j in range(30):
			t1 = time.time()
			sol = stochastic_greedy(i,15)[1]
			t2 = time.time()
			ls_times.append(t2-t1)
			ls_best_sol.append(sol)
			
			t1 = time.time()
			sol2 = stochastic_local_search(i,15)[1]
			t2 = time.time()
			gr_times.append(t2-t1)
			gr_best_sol.append(sol2)
		ls_avg_times.append(np.mean(ls_times))
		#np.std
		ls_avg_sol.append(np.mean(ls_best_sol))
		gr_avg_times.append(np.mean(gr_times))
		gr_avg_sol.append(np.mean(gr_best_sol))
		print(i)
	x = list(range(1,11))
	plt.subplot(121)
	plt.plot(x, ls_avg_times, 'c', label='LS')
	plt.plot(x, gr_avg_times, 'r', label="Greedy")
	plt.title("Avg times (on 30 runs) for len 15, number of words ranging ranging LS vs GR")
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	
	plt.subplot(122)
	plt.plot(x, ls_avg_sol, 'c', label="LS")
	plt.plot(x, gr_avg_sol, 'r', label="Greedy")
	plt.title("Avg constraints unfulfilled(on 30 runs) for len 15, number of words ranging LS vs GR")
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	
	plt.show()


# find_best_solutions()
