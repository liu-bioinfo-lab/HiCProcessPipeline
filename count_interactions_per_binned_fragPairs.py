#!/usr/bin/env python
'''
Created on Nov 30, 2011
revised on Jun 05, 2012

@author: ferhat
@edit : gunhan
'''
import sys
import math

USAGE = """USAGE:count_interactions_per_binned_fragPairs.py <pairsFileName> <REmidsFileName> <fragsize> <mappabilityThreshold> <distanceThreshold> <interactionsFileName> <mappableREmidsFilename>

Counts the number of interactions for each fragment pair, does the labeling of fragments (mappable(1)/not-mappable(0)) by using <mappabilityThreshold> value and outputs two files: <interactionsFileName> <mappableREmidsFilename>

<pairsFileName> is the main input file that contains interactions in this format:
"FCD059NACXX:5:2107:18137:185323 +       chr10   100000018       GCACCACCATACTCCAGCCTGGGCAACTGAGCAAGACTCCATCTCAAAA       +       chr5    66808096        ATCTTCTTGACTCCAGTCTTTCTTCTTATCAAAGACTGGATGTATAGGG"

<REmidsFileName> is the name of the input file that contains all the RE fragment midpoints as lines of the form 'chrno	REfragMidpoint'

<fragsize> size of the fragments that we bin the DNA segments into. Using 10000 for yeast, 100000 for human for now.

<mappabilityThreshold> n to take all RE fragments that have at least n valid interactions into account for downstream analysis.

<interactionsFileName> this file will have the output in the format 'chr1       fragmid1  chr2      fragmid2  noofinteractions'

<mappableREmidsFilename> this file is mainly for 'fit-hic' analysis and has some unnecessary columns. the format is:
'chrno       0       fragmid    number-of-interactions-this-fragment-participates-in       mappable-or-not(0/1)'

!! SORTING ACCORDING TO CHROMOSOMES and THEN MIDS IS DONE WTIH THIS SCRIPT (e.g. chr21 mid1 chr5 mid2 is converted to 5 mid2 21 mid1). 

"""

def countInteractionAndEliminateNonMappables(pairsFileName, REmidsFileName, fragsize, mappabilityThreshold, distanceThreshold, interactionsFileName, mappableREmidsFilename):
	#open all files
	REmidsFile=open(REmidsFileName, 'r')
	pairsFile =open(pairsFileName, 'r')
	outFileInteractions = open(interactionsFileName, 'w')
	outFileREmids = open(mappableREmidsFilename, 'w') 
	
	# to count the interactions between bins
	countInteractionsDict={}
	
	REmidsPerChr=[]
	chrList=[]
	#convert chromosome referencing to dictionary
	#allREmidsInteractionCount=[]
	#allREmids=[]

	#BEGIN
	#start reading all bin midpoints- assuming fixed bin sizes, which are each "fragsize" wide
	allREmids={}
	#line=REmidsFile.readline()
	#words=line.rstrip().split()
	#currChrName=words[0]
	currChrName="tmp_dummy_char"
	#midIndex=int(words[1])
	allREmidsInteractionCount={}
	for line in REmidsFile:
		words=line.rstrip().split()
		chrName=words[0]
		midIndex=int(float(words[1]))
		if currChrName!=chrName:
			if currChrName=="tmp_dummy_char":
				currChrName = chrName
				REmidsPerChr=[]
				#continue
			else:
				allREmids[currChrName]=REmidsPerChr
				allREmidsInteractionCount[currChrName]=[0 for _ in range(len(REmidsPerChr))]
				chrList.append(currChrName)
				print(currChrName +"\t"+str(len(REmidsPerChr))) 
				currChrName = chrName
				REmidsPerChr=[]
		REmidsPerChr.append(midIndex)
	#take care of the final chromosome here
	#allREmids.append(REmidsPerChr)
	allREmids[chrName]=REmidsPerChr
	#allREmidsInteractionCount.append([0 for _ in range(len(REmidsPerChr))])
	allREmidsInteractionCount[chrName]=[0 for _ in range(len(REmidsPerChr))]
	chrList.append(currChrName)
	print(currChrName +"\t"+str(len(REmidsPerChr))) 
	REmidsFile.close()
	#reading bins end here
	#END
	
	
	closePairs=0
	validPairs=0
	sameFragPairs=0
	
	#BEGIN
	#start reading the pairs file
	for inline in pairsFile:
		# one line:
		#FCD059NACXX:5:2107:18137:185323 + chr10 100000018 GCACCACCATACTCCAGCCTGGGCAACTGAGCAAGACTCCATCTCAAAA + chr5 66808096 ATCTTCTTGACTCCAGTCTTTCTTCTTATCAAAGACTGGATGTATAGGG
		tokenizedStr=inline.rstrip().split()

		readID=tokenizedStr[0]
		
		#convert chromosome referencing to dictionary
		#tempchr1=(tokenizedStr[2])[3:] # omit the 'chr' part
		#tempchr2=(tokenizedStr[6])[3:]# omit the 'chr' part
		chrName1=tokenizedStr[2]
		chrName2=tokenizedStr[6]

		#get rid of unwanted chromosomes, anything not found in mids file is unwanted
		if (chrName1 not in chrList) or (chrName2 not in chrList):
			continue

		#chrNo1= int(tempchr1)
		#chrNo2= int(tempchr2)
		#use chrName1 and chrName2 instead

		strandSign1= tokenizedStr[1]
		strandSign2= tokenizedStr[5]
		readStart1=int(tokenizedStr[3])
		readStart2 = int(tokenizedStr[7])
		REfragMidIndex1 = math.floor(readStart1 / fragsize)
		REfragMidIndex2 = math.floor(readStart2 / fragsize)
		#omit seq1
		#omit seq2


		#if chrNo1==chrNo2 and abs(readStart1-readStart2) < distanceThreshold:
		if chrName1==chrName2 and abs(readStart1-readStart2) < distanceThreshold:
			closePairs +=1
			continue

		#if chrNo1==chrNo2 and abs(REfragMidIndex1-REfragMidIndex2) <1: #eliminate interactions on same fragment
		if chrName1==chrName2 and abs(REfragMidIndex1-REfragMidIndex2) <1: #eliminate interactions on same fragment
			sameFragPairs +=1
			# continue
		
		#if not too close or same, then it is valid!
		validPairs +=1

		# count the number of interaction that each RE frag participates in
		allREmidsInteractionCount[chrName1][REfragMidIndex1] +=1
		#print(str(REfragMidIndex2)+"\t"+str(len(allREmidsInteractionCount[chrName2])))
		allREmidsInteractionCount[chrName2][REfragMidIndex2] +=1

		#if chrNo1<chrNo2 or (chrNo1==chrNo2 and REfragMidIndex1 < REfragMidIndex2):
		#	dictkey=str(chrNo1)+'.'+str(REfragMidIndex1)+'-'+str(chrNo2)+'.'+str(REfragMidIndex2)
		#else:
		#	dictkey= str(chrNo2)+'.'+str(REfragMidIndex2)+'-'+str(chrNo1)+'.'+str(REfragMidIndex1)
		if chrName1<chrName2 or (chrName1==chrName2 and REfragMidIndex1 < REfragMidIndex2):
			dictkey=chrName1+'.'+str(REfragMidIndex1)+'-'+chrName2+'.'+str(REfragMidIndex2)
		else:
			dictkey=chrName2+'.'+str(REfragMidIndex2)+'-'+chrName1+'.'+str(REfragMidIndex1)
		# dictkey will have two isoforms. If one of them exists we will increase the count, otherwise we will add dictkey1 to our dictionary

		# count the number of different isoforms (times) that interaction appears
		if dictkey in countInteractionsDict:
			oldval = countInteractionsDict.get(dictkey)
			countInteractionsDict[dictkey] = int(oldval)+1
		else: # neither exists
			countInteractionsDict[dictkey] = 1

	pairsFile.close()
	#all pairs read, closing file
	#END


	print(repr(closePairs) +"	interactions were discarded due to distance threshold of "+ repr(distanceThreshold))
	print(repr(sameFragPairs) +"	interactions were discarded since they happened on the same fragment but were further than "+ repr(distanceThreshold))
	print(repr(validPairs) +"	interactions were counted as valid")


	#print(str(chrList))
	#print(str(allREmidsInteractionCount.keys()))
	# print these numbers in a file that is ready for fit-hic analysis
	#for x in range(len(allREmidsInteractionCount)):
	for x in chrList:
		for y in range(len(allREmidsInteractionCount[x])):
			mappable=0
			if allREmidsInteractionCount[x][y] >= mappabilityThreshold:
				mappable=1
				
			outFileREmids.write(x+'\t0\t'+str(allREmids[x][y])+'\t'+str(allREmidsInteractionCount[x][y])+'\t'+str(mappable)+'\n') # ready for fit-hic analysis
	outFileREmids.close()


	for x in countInteractionsDict:
		tokenized = x.split('-')
		firstEnd = tokenized[0].split('.')
		secondEnd=  tokenized[1].split('.')
		chr1= firstEnd[0]
		mid1Index= int(firstEnd[1])
		chr2= secondEnd[0]
		mid2Index = int(secondEnd[1])
		#mid1 = allREmids[int(chr1)-1][mid1Index]
		#mid2 = allREmids[int(chr2)-1][mid2Index]
		mid1 = allREmids[chr1][mid1Index]
		mid2 = allREmids[chr2][mid2Index]

		#if allREmidsInteractionCount[int(chr1)-1][mid1Index] >= mappabilityThreshold and allREmidsInteractionCount[int(chr2)-1][mid2Index] >= mappabilityThreshold:
		#	outFileInteractions.write(chr1+'\t'+str(mid1)+'\t'+chr2+'\t'+str(mid2)+'\t'+str(countInteractionsDict.get(x))+'\n') # ready for Fig2a and Heatmaps
		if allREmidsInteractionCount[chr1][mid1Index] >= mappabilityThreshold and allREmidsInteractionCount[chr2][mid2Index] >= mappabilityThreshold:
			outFileInteractions.write(chr1+'\t'+str(mid1)+'\t'+chr2+'\t'+str(mid2)+'\t'+str(countInteractionsDict.get(x))+'\n') # ready for Fig2a and Heatmaps

	outFileInteractions.close()
	
	
if __name__ == "__main__":
	if (len(sys.argv) != 8):
	  sys.stderr.write(USAGE)
	  sys.exit(1)
	countInteractionAndEliminateNonMappables(str(sys.argv[1]),str(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),str(sys.argv[6]),str(sys.argv[7]))
	
