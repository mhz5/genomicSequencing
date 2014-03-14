#!/usr/bin/env python

import time

# Bar codes for identifying relevant lines
chr4bar1	= "GAGCGA"
chr4bar2	= "GATGTG"
CFTRbar1	= "CGGAAT"
CFTRbar2	= "CGTGGC"
MCFTRbar1	= "TGCGTA"
MCFTRbar2	= "TTCTGG"

# ==== Begin sequence constants ====

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
interestMouse = "AGAAAATATCATTGG"     # Portion of interest (same as CFTR)

origMouse  = "TCTGCTCTCAATTTTCTTGGATTATGCCGGGTACTATCAA" 
origMouse += interestMouse 
origMouse += "TGTTTCCTATGATGAGTACAGATATAAGAGTGTTGTCAAAGCTTGCC"

# Desired product post-mutation.
mutInterestMouse = "GGAGAACATTATCTTTGG"    # Mutation portion of interest (same as CFTR)

mutMouse  = "TCTGCTCTCAATTTTCTTGGATTATGCCGGGTACTATCAA" 
mutMouse += mutInterestMouse 
mutMouse += "TGTTTCCTATGATGAGTACAGATATAAGAGTGTTGTCAAAGCTTGCC"

# Chr4
interestChr4 = "ATTAAAGAAAATATCCTATTTGGTCATTTATCATGAAGATA"
interestChr4Trunc = "TTTGGTCATTTATCATGAAGATA"

origChr4  = "TTCTACTAAAAGAAAACTTCTGTGTCCTA" 
origChr4 += interestChr4 
origChr4 += "ATGATAAATGTTAAATAACGTCTTGTTTGCATTAAGTCTGTGGGATG"

# ==== End sequence constants ====

infile = open("Sample_021714-1_005/021714-1_005_ACAGTG_L001_R1_003.fastq", "r")

startTime = time.time()

count = 0
chr4bar1	= "GAGCGA"
chr4bar2	= "GATGTG"
CFTRbar1	= "CGGAAT"
CFTRbar2	= "CGTGGC"
MCFTRbar1	= "TGCGTA"
MCFTRbar2

for line in infile:
	bar = line[0:6]
	if (bar == chr4bar1 or
		bar == chr4bar2 or
		bar == CFTRbar1 or
		bar == CFTRbar2 or
		bar == MCFTRbar1 or
		bar == MCFTRbar2):

		print count
		print line
		count += 1
	
print "\n\n\n======================"
print "Total time:"
print time.time() - startTime