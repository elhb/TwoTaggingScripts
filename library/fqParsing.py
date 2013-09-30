######################### MAIN #########################

def lib_main():
	pass

######################### FUNCTIONS #########################

def comp(str):
	return str.replace("A","X").replace("T","A").replace("X","T").replace("G","X").replace("C","G").replace("X","C")

def revcomp(str):
	return comp(str[::-1])

def getPairs(r1,r2):
	""" yield one readpair at a time from indata
	"""
	import re

	import sys

	# set the counters to initial values
	counter = 0
	tmp=0
	header="NA"

	file1 = r1
	file2 = r2
	
	#check if files are gzipped
	if file1.split('.')[-1] in ['gz','gzip']:
		import gzip
		file1 = gzip.open(file1)
	else: file1 = open(file1,'r')
	if file2.split('.')[-1] in ['gz','gzip']:
		import gzip
		file2 = gzip.open(file2)
	else: file2 = open(file2,'r')

	# itarate through fastq files and return readpairs
	for r1line in file1:
		r2line = file2.readline() #get line from  read2 file
		
		tmp+=1 # increment linecount
		
		# depending on line number (within entry) do ...	
		if tmp == 1: #header check match between files
			counter+=1 # increase entry counter
			header=r1line
			if r1line.split(" ")[0] != r2line.split(" ")[0]:
				sys.stderr.write('Error mismatching headers!');
				raise ValueError
		elif tmp == 2: # sequence check that its DNA and save sequences till later
			if counter == 1:
				sys.stderr.write('Checking data type of read 1 in pair 1 ... ')
				match = re.match("^[AGTCN]+$",r1line.rstrip())
				if match: sys.stderr.write('this is DNA data.\n')
				else:
					sys.stderr.write(' this is not a DNA sequence ('+r1line.rstrip()+') could something be wrong with your fastq file?.\n\n');
					raise ValueError
			r1seq = r1line.rstrip()
			r2seq = r2line.rstrip()
		elif tmp == 3: # "+"-line do some format check
				if counter in {1:True,67:True,438:True,9675:True,53678:True,864513:True,1337354:True,317955:True,1226844:True,20389:True,118261:True}:
					if r1line[0] != r2line[0] or r1line[0] != '+': sys.stderr.write('Error Format not fastq!');raise ValueError#os.kill(MASTER);sys.exit(1);#REALLYNOTOPTIMAL
		elif tmp == 4: # quality values and end of entry, reset counter and yeild readpair
				tmp=0 # reset line counter
				r1qual = r1line.rstrip() #get qual strings
				r2qual = r2line.rstrip()
				
				#yield readpair
				yield readpair(header.rstrip(), read(header.rstrip(),r1seq,r1qual), read(header.rstrip(),r2seq,r2qual))

def hamming_distance(s1, s2):
	assert len(s1) == len(s2), 'Error: '+str(len(s1)) + ' != ' + str(len(s2))
	#return sum(ch1 not in uipac(ch2) for ch1, ch2 in zip(s1, s2))
	mms = 0
	for ch1, ch2 in zip(s1, s2):
		if ch1 == 'N' or ch2 == 'N': mms+=1;continue # N will always count as missmatch
		match = False
		for b1 in uipac(ch1, back='bases'):
			for b2 in uipac(ch2, back='bases'):
				if b1 == b2: match = True
		if not match: mms +=1
	return mms

def bufcount(filename):
	""" returns the number of lines in a file
	"""
	import gzip
	if filename.split('.')[-1] in ['gz','gzip']: f = gzip.open(filename)
	else: f = open(filename)
	lines = 0
	buf_size = 1024 * 1024
	read_f = f.read # loop optimization
	
	buf = read_f(buf_size)
	while buf:
		lines += buf.count('\n')
		buf = read_f(buf_size)
		f.close
	return lines

def uipac(bases, back='uipac'): #U	Uracil NOT SUPPORTED!!!
	if back == 'uipac':
		if 'N' in bases: return 'N'
		uniqbases={}
		for i in bases:
			uniqbases[i]=True
		bases = uniqbases.keys()
		if 'U' in bases: sys.stderr.write('WARNING in function "uipac(bases)": Uracil NOT SUPPORTED!')
		if len(bases)==1:
			if 'A' in bases: return 'A' #A	Adenine
			if 'C' in bases: return 'C' #C	Cytosine
			if 'G' in bases: return 'G' #G	Guanine
			if 'T' in bases: return 'T' #T	Thymine
			#U	Uracil NOT SUPPORTED!!!
		if len(bases)==2:
			if 'A' in bases and 'G' in bases: return 'R' #R	Purine (A or G)
			if 'C' in bases and 'T' in bases: return 'Y' #Y	Pyrimidine (C, T, or U)
			if 'A' in bases and 'C' in bases: return 'M' #M	C or A
			if 'T' in bases and 'G' in bases: return 'K' #K	T, U, or G
			if 'A' in bases and 'T' in bases: return 'W' #W	T, U, or A
			if 'C' in bases and 'G' in bases: return 'S' #S	C or G
		if len(bases)==3:
			if 'C' in bases and 'T' in bases and 'G' in bases: return 'B' #B	C, T, U, or G (not A)
			if 'A' in bases and 'T' in bases and 'G' in bases: return 'D' #D	A, T, U, or G (not C)
			if 'A' in bases and 'T' in bases and 'C' in bases: return 'H' #H	A, T, U, or C (not G)
			if 'A' in bases and 'C' in bases and 'G' in bases: return 'V' #V	A, C, or G (not T, not U)
		if len(bases)==4:
			return 'N' #N	Any base (A, C, G, T, or U)
	elif back == 'bases':
		if   bases == 'R': return ['A','G'] 
		elif bases == 'Y': return ['C','T']
		elif bases == 'M': return ['A','C']
		elif bases == 'K': return ['G','T']
		elif bases == 'W': return ['A','T']
		elif bases == 'S': return ['C','G']
		elif bases == 'B': return ['C','T','G']
		elif bases == 'D': return ['A','T','G']
		elif bases == 'V': return ['A','C','G']
		elif bases == 'H': return ['A','C','T']
		elif bases == 'N': return ['A','G','T','C']
		elif 'A' == bases: return 'A' #A	Adenine
		elif 'C' == bases: return 'C' #C	Cytosine
		elif 'G' == bases: return 'G' #G	Guanine
		elif 'T' == bases: return 'T' #T	Thymine


def UIPAC2REGEXP(string):
    return string.replace('R','[AG]').replace('Y','[CT]').replace('S','[GC]').replace('W','[AT]').replace('K','[GT]').replace('M','[AC]').replace('B','[CGT]').replace('D','[AGT]').replace('H','[ACT]').replace('V','[ACG]').replace('N','.')

######################### CLASSES #########################

class sequence():
	def __init__(self,header,seq,qual):
		self.header = header.rstrip()
		self.qual = qual.rstrip()
		self.seq = seq.rstrip()
		assert len(self.qual) == len(self.seq), 'Error: qual and seq has different lengths!\n'
		self.len = len(seq.rstrip())
	
	def subseq(self,start,end):
		return sequence(self.header,self.seq[start:end],self.qual[start:end])
	
	def revcomp(self):
		''' Takes a sequence and reversecomplements it'''
		complementary = self.comp()
		return sequence(complementary.header,complementary.seq[::-1],complementary.qual[::-1])
	
	def comp(self):
		''' Takes a sequence and complements it'''
		complement = {'A':'T','T':'A',
					  'C':'G','G':'C',
					  'N':'N',
					  'R':'Y','Y':'R',
					  'K':'M','M':'K',
					  'B':'V','V':'B',
					  'D':'H','H':'D',
					  }
		compseq = "".join([complement.get(nt.upper(), '') for nt in self.seq])
		return sequence(self.header,compseq,self.qual)

class read(sequence):
    "Represents one of several reads from a DNA fragment"

class readpair():
	""" object representing an illumina cluster """

	def __init__(self,header,r1,r2):
		self.header = header
		self.r1 = r1 #first read
		self.r2 = r2 #second read
		self.fwdPrimer = sequence('fwdPrimer','GATGGAMTCCGRTTTAKYGGCAT','GATGGAMTCCGRTTTAKYGGCAT')
		self.revPrimer = sequence('revPrimer','ATCTGACATAAAGATCTTGACC','ATCTGACATAAAGATCTTGACC')

	def identify(self, handle, config):
		[handle_start, handle_end] = self.matchHandle(handle, config, self.r1)
		self.handle_start = handle_start
		self.handle_end   = handle_end
		return 0

	def identifyIllumina(self, config):
		#handle = sequence('illuminaUniversal','AGATCGGAAGAGC','AGATCGGAAGAGC')
		config.chandlemissmatch = 3
		#handle = sequence('illumina','AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC','AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC')
		handle = sequence('illumina','AGATCGGAAGAGCACACGTCT','AGATCGGAAGAGCACACGTCT')
		[handle_start, handle_end] = self.matchHandle(handle, config, self.r1)
		if handle_start: self.r1.illuminaadapter = True
		else: self.r1.illuminaadapter = False
		
		#handle = sequence('illumina','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
		handle = sequence('illumina','AGATCGGAAGAGCGTCGTGT','AGATCGGAAGAGCGTCGTGT')
		[handle_start, handle_end] = self.matchHandle(handle, config, self.r2)
		if handle_start: self.r2.illuminaadapter = True
		else: self.r2.illuminaadapter = False
		
		if self.r1.illuminaadapter or self.r2.illuminaadapter:
			self.isillumina = True
		else:	self.isillumina = False
		return 0

	def matchPrimers(self, matchfunk=hamming_distance, primermissmatch=0,tagmissmatch=0):

		self.fwdPrimer.read	= None
		self.fwdPrimer.start	= None
		self.fwdPrimer.end	= None
		self.revPrimer.read	= None
		self.revPrimer.start	= None
		self.revPrimer.end	= None

		import re
		for primer in [self.fwdPrimer, self.revPrimer]:
			
			for read in [self.r1, self.r2]:
				
				perfect_match = re.search(UIPAC2REGEXP(primer.seq), read.seq)
				if perfect_match:
				    primer.start = perfect_match.start()
				    primer.end = perfect_match.end()
				    primer.read = read
				    primer.missmatch = 0
				    break
				elif primermissmatch:
					mindist = [10000,-1]
					for i in range(len(read.seq)):
						if i < 8: continue
						if len(read.seq)-i > len(primer.seq):
							dist = matchfunk(primer.seq,read.seq[i:i+len(primer.seq)])
						else: dist = 1000
						if dist < mindist[0]: mindist =[dist,i]
						if i > 8: break
					if mindist[0] <= primermissmatch:
						i = mindist[1]
						primer.missmatch = mindist[0]
						primer.start = i
						primer.end = i+len(primer.seq)
						primer.read = read
		
		if self.fwdPrimer.read:
			self.fwdTagSeq = self.fwdPrimer.read.seq[:self.fwdPrimer.start]
			self.fwdWell = None
			if len(self.fwdTagSeq) > 8: self.fwdTagSeq = '>8bp'
			elif len(self.fwdTagSeq)<8: self.fwdTagSeq = '<8bp'
			else:
				try: self.fwdWell = tags[self.fwdTagSeq]
				except KeyError:
					for tag in tags:
						dist = hamming_distance(tag,self.fwdTagSeq)
						if dist <= tagmissmatch:
							self.fwdTagmissmatch = dist
							self.fwdWell = tags[tag]
							break
		if self.revPrimer.read:
			self.revTagSeq = self.revPrimer.read.seq[:self.revPrimer.start]
			self.revWell = None
			if len(self.revTagSeq) > 8: self.revTagSeq = '>8bp'
			elif len(self.revTagSeq)<8: self.revTagSeq = '<8bp'
			else:
				try: self.revWell = tags[self.revTagSeq]
				except KeyError:
					for tag in tags:
						dist = hamming_distance(tag,self.revTagSeq)
						if dist <= tagmissmatch:
							self.revTagmissmatch = dist
							self.revWell = tags[tag]
							break
		if self.fwdPrimer.read and self.revPrimer.read: self.effectivelength = len(self.revPrimer.read.seq[self.revPrimer.end:])+len(self.fwdPrimer.read.seq[self.fwdPrimer.end:])
		else:self.effectivelength = 0
		return 0

tags = {
    'TCTCTGTG':'A1',
    'TGTACGTG':'A2',
    'ATCGTCTG':'A3',
    'TAGCTCTG':'A4',
    'AGTATCTG':'A5',
    'TCGAGCTG':'A6',
    'TCATACTG':'A7',
    'TACGACTG':'A8',
    'ACTCACTG':'A9',
    'AGAGTATG':'A10',
    'AGCTGATG':'A11',
    'TATCGATG':'A12',
    'ATGCGATG':'B1',
    'ACGTCATG':'B2',
    'TCATGTCG':'B3',
    'TAGCGTCG':'B4',
    'TCTACTCG':'B5',
    'ATGACTCG':'B6',
    'ATCTATCG':'B7',
    'ACAGATCG':'B8',
    'ATACTGCG':'B9',
    'TATATGCG':'B10',
    'TGCTCGCG':'B11',
    'ATCGCGCG':'B12',
    'TAGTAGCG':'C1',
    'AGATAGCG':'C2',
    'TGTGAGCG':'C3',
    'TCACAGCG':'C4',
    'ACTGTACG':'C5',
    'TGCGTACG':'C6',
    'TCGCTACG':'C7',
    'TACTGACG':'C8',
    'AGACGACG':'C9',
    'TGTAGACG':'C10',
    'ACGAGACG':'C11',
    'ATATCACG':'C12',
    'TCAGCACG':'D1',
    'TAGACACG':'D2',
    'AGCACACG':'D3',
    'ATGTGTAG':'D4',
    'ACTCGTAG':'D5',
    'TGCAGTAG':'D6',
    'TGATCTAG':'D7',
    'TACGCTAG':'D8',
    'TCGTATAG':'D9',
    'AGACATAG':'D10',
    'AGCGTGAG':'D11',
    'ATGATGAG':'D12',
    'ACATCGAG':'E1',
    'TCTGCGAG':'E2',
    'ATAGAGAG':'E3',
    'TATCAGAG':'E4',
    'ACGCAGAG':'E5',
    'ACAGTCAG':'E6',
    'TCTATCAG':'E7',
    'TAGTGCAG':'E8',
    'TGACGCAG':'E9',
    'ATCAGCAG':'E10',
    'TGCTACAG':'E11',
    'AGTGACAG':'E12',
    'ACTGTGTC':'F1',
    'TACATGTC':'F2',
    'ATGACGTC':'F3',
    'AGCGAGTC':'F4',
    'TCGCAGTC':'F5',
    'ATACAGTC':'F6',
    'TGCGTCTC':'F7',
    'TCACTCTC':'F8',
    'ATCTGCTC':'F9',
    'TGTAGCTC':'F10',
    'ACGTACTC':'F11',
    'TCTGACTC':'F12',
    'ACGCTATC':'G1',
    'ATCATATC':'G2',
    'TCGTGATC':'G3',
    'TGACGATC':'G4',
    'TGCTCATC':'G5',
    'TATGCATC':'G6',
    'ACAGCATC':'G7',
    'AGTACATC':'G8',
    'AGTGCTGC':'G9',
    'TGCGATGC':'G10',
    'ATGCATGC':'G11',
    'TCACATGC':'G12',
    'AGAGTCGC':'H1',
    'ACTATCGC':'H2',
    'TAGATCGC':'H3',
    'TCATGCGC':'H4',
    'TACTACGC':'H5',
    'ATATACGC':'H6',
    'TGTCACGC':'H7',
    'AGTCTAGC':'H8',
    'ATGTGAGC':'H9',
    'TAGCGAGC':'H10',
    'ACACGAGC':'H11',
    'TCTAGAGC':'H12'
    }

if __name__ == "__main__": lib_main()

