#!/usr/bin/env python

import time
from copy import copy

# Bar codes for identifying relevant lines
chr4bar1	= "GAGCGA"
chr4bar2	= "GATGTG"
CFTRbar1	= "CGGAAT"
CFTRbar2	= "CGTGGC"
MCFTRbar1	= "TGCGTA"
MCFTRbar2	= "TTCTGG"

# === Begin sequence constants ===

# CFTR (what does WT stand for??)
interestCFTR = "AGAAAATATCATTGG"      # Portion of interest (same as mouse)

origCFTR  = "TGTTCTCAGTTTTCCTGGATTATGCCTGGCACCATTAA"
origCFTR += interestCFTR
origCFTR += "TGTTTCCTATGATGAATATAGATACAGAAGCGTCATCAAAGCATGCCA"

# Desired product post-mutation
mutInterestCFTR = "GGAGAACATTATCTTTGG"      # Mutation portion of interest (same as mouse)

mutCFTR  = "TGTTCTCAGTTTTCCTGGATTATGCCTGGCACCATTAA"
mutCFTR += mutInterestCFTR
mutCFTR += "TGTTTCCTATGATGAATATAGATACAGAAGCGTCATCAAAGCATGCCA"

# Mouse
interestMCFTR = "AGAAAATATCATTGG"     # Portion of interest (same as CFTR)
# Take 11 char slice.
origMCFTR  = "TCTGCTCTCAATTTTCTTGGATTATGCCGGGTACTATCAA" 
origMCFTR += interestMCFTR 
origMCFTR += "TGTTTCCTATGATGAGTACAGATATAAGAGTGTTGTCAAAGCTTGCC"

# Desired product post-mutation.
mutInterestMCFTR = "GGAGAACATTATCTTTGG"    # Mutation portion of interest (same as CFTR)

mutMCFTR  = "TCTGCTCTCAATTTTCTTGGATTATGCCGGGTACTATCAA" 
mutMCFTR += mutInterestMCFTR 
mutMCFTR += "TGTTTCCTATGATGAGTACAGATATAAGAGTGTTGTCAAAGCTTGCC"

# CHR4
interestCHR4 = "ATTAAAGAAAATATCCTATTTGGTCATTTATCATGAAGATA"
interestCHR4Trunc = "TTTGGTCATTTATCATGAAGATA"

origCHR4  = "TTCTACTAAAAGAAAACTTCTGTGTCCTA" 
origCHR4 += interestCHR4 
origCHR4 += "ATGATAAATGTTAAATAACGTCTTGTTTGCATTAAGTCTGTGGGATG"

# === End sequence constants ===

# === Begin output variables ===

# Enum to hold experiment types.
class E:
	CFTR 	= 0
	MCFTR 	= 1
	CHR4 	= 2

statList = {'unmutROI': "",
			'mutInterest': "",	# Portion to be replaced.
			'origSeq': "",		# Original sequence (to compare)
			'mutSeq': "",		# Expected mutated sequence (to compare)
			'total': 0,			# Total num of sequences from specified experiment.
			'unmut': 0,			# Unmutated sequences.
			'mut': 0,			# Mutated sequences.
			'noMatch': 0,		# No match sequences.
			'mutOrig': 0,		# Num mutations in original portion of interest.
			'flankOrig': 0,		# Num mutations in original flanking regions.
			'mutMut': 0,		# Num mutations in MUTATED portion of interest.
			'flankMut': 0 		# Num mutations in MUTATED flanking regions.
			}

# Object holding analysis variables.
var = {E.CFTR: copy(statList), E.MCFTR: copy(statList), E.CHR4: copy(statList)}

# Store known data
var[E.CFTR]['unmutInterest'] = interestCFTR
var[E.MCFTR]['unmutInterest'] = interestCFTR
var[E.CHR4]['unmutInterest'] = interestCFTR

var[E.CFTR]['mutInterest'] 	= mutInterestCFTR
var[E.MCFTR]['mutInterest'] = mutInterestMCFTR
var[E.CHR4]['mutInterest'] 	= interestCHR4

var[E.CFTR]['origSeq']		= origCFTR
var[E.MCFTR]['origSeq']		= origMCFTR
var[E.CHR4]['origSeq']		= origCHR4

var[E.CFTR]['mutSeq']		= mutCFTR
var[E.MCFTR]['mutSeq']		= mutMCFTR
var[E.CHR4]['mutSeq']		= origCHR4

# === End output variables ===

# === Utility functions ===

# This is a huge drain.
def strMatches(a, b):
	if a == b:
		return len(a)
	minLen = min(len(a), len(b))
	res 
	for i in range(0, minLen):
		if a[i] == b[i]:
			res += 1
	return res

# Return > match score of two sequences of same length
# If not same length, only process up to the length of the shorter sequence.
def matchScore(seq1, seq2):
	match = 1 		# match
	mismatch = -1 	# mismatch
	gap = -1 	# gap

	n = min(len(seq1), len(seq2))
	# Initialize dynamic programming array
	arr = [[0 for j in range(n + 1)] for i in range(n + 1)]

	maxScore = 0

	# Fill in array
	for i in range(1, n + 1):
		# Only consider a slice of width 6 to save time.
		for j in range(max(1, i - 3), min(n + 1, i + 3)):
			
			diag = arr[i - 1][j - 1]	# match/mismatch case

			if seq1[i - 1] == seq2[j - 1]:	# match
				diag += match
			else:					# mismatch
				diag += mismatch

			gap1 = arr[i - 1][j] + gap 	# deletion case
			gap2 = arr[i][j - 1] + gap 	# insertion case

			arr[i][j] = max(0, diag, gap1, gap2)

			maxScore = max(maxScore, arr[i][j])
	
	return maxScore

# Return => void
# Primary processing function
def process(line, exp, start, end, scoreExpected):
	var[exp]['total'] += 1
	interest = line[start:end]
	
	scoreTotal = 0
	scoreROI = matchScore(interest, var[exp]['mutInterest'])
	scoreDiff = scoreExpected - scoreROI

	if scoreDiff < 8:			# Mutated
		print "cool region: " + interest
		print "actual     : " + var[exp]['mutInterest']
		print scoreROI
		print exp
		var[exp]['mut'] += 1
		var[exp]['mutMut'] += scoreDiff
		scoreTotal = matchScore(line, var[exp]['mutSeq'][0:len(line)])
		var[exp]['flankMut'] += len(line) - scoreTotal - var[exp]['mutMut']

	else:								# Unmutated
		var[exp]['unmut'] += 1
		scoreROI = matchScore(interest, var[exp]['unmutInterest'])
		var[exp]['mutOrig'] += len(var[exp]['unmutInterest']) - scoreROI
		scoreTotal = matchScore(line, var[exp]['origSeq'][0:len(line)])
		var[exp]['flankOrig'] += len(line) - scoreTotal - var[exp]['mutOrig']



# === End Utility functions ===

infile = open("testInput", "r")

startTime = time.time()

# Identify experiment by barcode.
chr4bar1	= "GAGCGA"
chr4bar2	= "GATGTG"

CFTRbar1	= "CGGAAT"
CFTRbar2	= "CGTGGC"

MCFTRbar1	= "TGCGTA"
MCFTRbar2	= "TTCTGG"

mice = 0
mutMice = 0
unmutMice = 0

for line in infile:
	bar = line[0:6]
	# Chr4
	# (strMatches(bar, chr4bar1) > 4 or
	# 	strMatches(bar, chr4bar2) > 4):

# NOTE: Give more room for error in the ROI
	if (bar == chr4bar1 or
		bar == chr4bar2):
		process(line, E.CHR4, 29 + 6, 76, 41)

	# CFTR
	if (bar == CFTRbar1 or
		bar == CFTRbar2):
		process(line, E.CFTR, 38 + 6, 62, 18)	# errors
	
	# Mouse
	if (bar == MCFTRbar1 or
		bar == MCFTRbar2):
		process(line, E.MCFTR, 40 + 6, 64, 18)

# Output

print "\n\n\n======================"

print "CFTR total: " + str(var[E.CFTR]['total'])
print "CFTR mutated: " + str(var[E.CFTR]['mut'])
print "CFTR no-match: " + str(var[E.CFTR]['noMatch'])
print "CFTR ROI mutations (original): " + str(var[E.CFTR]['mutOrig'] / 2)
print "CFTR flanking mutations (original): " + str(var[E.CFTR]['flankOrig'] / 2)
print "CFTR ROI mutations (mutated): " + str(var[E.CFTR]['mutMut'] / 2)
print "CFTR flanking mutations (mutated): " + str(var[E.CFTR]['flankMut'] / 2)

print "MCFTR total: " + str(var[E.MCFTR]['total'])
print "MCFTR mutated: " + str(var[E.MCFTR]['mut'])
print "MCFTR no-match: " + str(var[E.MCFTR]['noMatch'])
print "MCFTR ROI mutations (original): " + str(var[E.MCFTR]['mutOrig'] / 2)
print "MCFTR flanking mutations (original): " + str(var[E.MCFTR]['flankOrig'] / 2)
print "MCFTR ROI mutations (mutated): " + str(var[E.MCFTR]['mutMut'] / 2)
print "MCFTR flanking mutations (mutated): " + str(var[E.MCFTR]['flankMut'] / 2)

print "Chr4 total: " + str(var[E.CHR4]['total'])
print "CHR4 mutated: " + str(var[E.CHR4]['mut'])
print "Chr4 ROI mutations: " + str(var[E.CHR4]['mutOrig'] / 2)
print "Chr4 flanking mutations: " + str(var[E.CHR4]['flankOrig'] / 2)
print "Chr4 no-match: " + str(var[E.CHR4]['noMatch'] / 2)

print "Total time: "
print time.time() - startTime