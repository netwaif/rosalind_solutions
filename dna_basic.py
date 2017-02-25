import re


def get_input(filename="input.txt"):
	f = file("input.txt", "r")
	input = ""
	for line in f:
		input += line.replace("\n", " ").replace("\r", "").replace("  ", " ")
	args = ()
	for a in input.split(" "):
		args += (a,)
	f.close()
	return args


def num_args(raw_args, args_count):
	args = []
	i = 1
	for a in raw_args:
		if not a.isdigit():
			continue
		elif i <= args_count:
			args.append(int(a))
			++i
		else:
			break
	return args


def dna_stat(dna_str, rna=False):
	"""
	:param dna_str: string of dna/rna
	:param rna: True for RNA
	:return: dict of numbers for each ACGT + key "sum" for total # of letters
	"""
	stat = {}
	rna = dna_str.find("U") or rna
	if rna == False:
		stat = dict(zip((n for n in "ACGT"), (dna_str.count(n) for n in "ACGT")))
		stat["rna"] = False
	elif rna == True:
		stat = dict(zip((n for n in "ACGU"), (dna_str.count(n) for n in "ACGU")))
		stat["rna"] = True
	stat["sum"] = sum(i for i in stat.values())
	return stat


def dna_reverse_complement(dna_str):
	res_seq = ""
	matrix = {"A": "T", "T": "A", "C": "G", "G": "C"}
	for n in dna_str[::-1]:
		res_seq += matrix[n]
	return res_seq


def fib(generations, start=1, fertility=1, mortality=0):
	a, b = start, 0
	cache = {0: start, 1: start}
	for i in xrange(1, generations + 1):
		if i >= mortality:
			a, b = a + fertility * b - cache[i - mortality], a
			cache.__delitem__(i - mortality)
		else:
			a, b = a + fertility * b, a
		if b < 0:
			b = 0
		cache[i] = b
	cache.clear()
	return b


def dna_gc(dna_str):
	stat = dna_stat(dna_str)
	return 100.0 * (stat['C'] + stat['G']) / stat['sum']


def dna_read_FASTA(dna_hash={"local": True}, filename="input.txt"):
	re_ID = re.compile(r'\s*>(Rosalind_[0-9]{4})')
	re_DNA = re.compile(r'([^ACGT]{1})')

	with open(filename, "r") as f:
		current_ID = ''
		for line in f:
			if re_ID.match(line):
				m_ID = re_ID.match(line)
				current_ID = m_ID.group(1)
				dna_hash[current_ID] = ""
			elif current_ID != '':
				dna_hash[current_ID] += re_DNA.sub("", line.strip().upper())

		if dna_hash.has_key("local"):
			dna_hash.pop("local")
			return dna_hash


def dna_gc_base(dna_base, freq_base):
	curr_max = 0
	curr_max_ID = ""
	for id, dna in dna_base.iteritems():
		freq_base[id] = dna_gc(dna)
		if freq_base[id] > curr_max:
			curr_max = freq_base[id]
			curr_max_ID = id
	return curr_max_ID
