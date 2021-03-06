from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

#SARS-2 (Covid-19): NC_045512

# Contains 59 DNA strands from the coronavirus family
arr = [""]*59
strippedArr = [""]*59
# For array index
i = 0

# For creating sars-2 primers
for seq_record in SeqIO.parse("seqAlign-59.fasta", "fasta"):
    arr[i] = str(seq_record.seq).lower()
    strippedArr = arr[i].replace("-","")
    i += 1

# For interested sequence
ignore = 23
intSeq = str(arr[ignore])
# stripped version
strippedSeq = intSeq.replace("-","")

# Candidate primers and locations
# Candidate primer set generated by initial notion of "unique"
#         and the set of 9 aligned sequences.
# candidatePrimers = {
#                     ('acatccactccctgtaatc', 'gcctaggatacaacacgtc'), 
#                     ('cgttatgctagtgctgttg', 'gcttggctttctgaacttg'), 
#                     ('aaatggctgtcgcttatgc', 'gcaatggccaatccatttc'), 
#                     ('gttactaggctgcatgatg', 'tgtacgtggttgtactacc'), 
#                     ('ttctggcatggacaccgcattg', 'cctgatggttgctgagagg')
#                    }
# candidatePrimerLocs = {
#                        (30339, 30473), 
#                        (23587, 24075), 
#                        (58, 186), 
#                        (2300, 2569), 
#                        (21488, 21894), 
#                        (29029, 29244)
#                       }


# Primers generated by new notion of "unique" and set of 59 sequences
candidatePrimers = {('gtgaggacaatcagacaac', 'cctctttaggaatctcagc')}
candidatePrimerLocs = {(3244, 4041)}
                      
# REQUIRES:
#   - nucStr is a valid DNA strand
# ENSURES:
#   - complement strand is designed
#       (3',5' end is unknown)
def genComplement(nucStr):
    build = ""
    for nucChar in nucStr:
        if (nucChar == 'a'):
            build += 't'
        elif (nucChar == 't'):
            build += 'a'
        elif (nucChar == 'c'):
            build += 'g'
        elif (nucChar == 'g'):
            build += 'c'
        else:
            raise Exception (nucChar + " is not valid DNA strand")
    return build

# For calculating TM score
def getTmScore(seq):
    score = 0
    for chr in seq:
        if chr == 'g' or chr == 'c':
            score += 4
        elif chr == 'a' or chr == 't':
            score += 2
        else:
            raise Exception (chr + " not valid DNA strand")
    return score

# Make sure they only appear in the area we want
def filterSeq(s, seq):
    for i in range(0, ignore):
        if str(strippedArr[i]).find(seq) != -1:
            return False
    for i in range(ignore + 1, len(arr)):
        if str(strippedArr[i]).find(seq) != -1:
            return False
    return s == strippedSeq.rfind(seq)

# Find location where primer most likely came from
def findLoc (seq, isRev):
    if isRev:
        seq = genComplement(seq)
        seq = seq[::-1]
    maxIndex = 0
    maxScore = -1
    for i in range(len(strippedSeq)-len(seq)):
        # Calculate score
        tempScore = 0
        for j in range(len(seq)):
            if seq[j] == strippedSeq[i+j]:
                # Extreme bias in weighting for matches at 3' end, since last
                #   5 nucleotides at 3' end drastically affect how well primer
                #       anneals to strand
                if isRev == False and j > len(seq) - 5:
                    tempScore += 1 # Bias for forward primer
                if isRev == True and j < 5:
                    tempScore += 1 # Bias for reverse primer
                tempScore += 1
        if tempScore > maxScore:
            maxScore = tempScore
            maxIndex = i
    return maxIndex

referencePrimers = set()

# Obtained from PrimerBLAST tool (input reference sequence from NCBI database)
referencePrimers.add(("CGGATGGCTTATTGTTGGCG", "TTGTGCTTACAAAGGCACGC"))
referencePrimers.add(("GCCGCTGTTGATGCACTATG", "CTCCAAGCAGGGTTACGTGT"))
referencePrimers.add(("TGTCGTTGACAGGACACGAG", "TTACCTTTCGGTCACACCCG"))
referencePrimers.add(("TGAGCCAGTGCTCAAAGGAG", "CGCCAACAATAAGCCATCCG"))
referencePrimers.add(("TACGGCGCCGATCTAAAGTC", "GCGACCACCCTTACGAAGAA"))
referencePrimers.add(("ATCTGTGTGGCTGTCACTCG", "ACTCGTGTCCTGTCAACGAC"))
referencePrimers.add(("AGAATGGAGAACGCAGTGGG", "TGAGAGCGGTGAACCAAGAC"))
referencePrimers.add(("TGTTCCTCGGAACTTGTCGG", "CTGAGGTGTGTAGGTGCCTG"))
referencePrimers.add(("GACTAATTCTCCTCGGCGGG", "TGTACCCGCTAACAGTGCAG"))
referencePrimers.add(("TCTTGGTTCACCGCTCTCAC", "TTGTTAGCAGGATTGCGGGT"))

# Tests for the reference primers
# Modified percentCG for string input
def pCG(x):
    ret = 0
    for i in x:
        if i == 'g' or i == 'c':
            ret += 1
    return ret/len(x)*100

# Modified findRepeats for string input
def fRepeats(x):
    s = set()
    s.add("aaaa")
    s.add("tttt")
    s.add("cccc")
    s.add("gggg")
    for i in range(len(x)-4):
        if intSeq[i:(i+4)] in s:
            return i
    return -1

# Kept calculating ranges, helper function for that
def inRange(x,lower,upper):
    return x >= lower and x <= upper

for (a,b) in referencePrimers:
    forward = a.lower()
    backward = b.lower()
    # Test to see if their primers meet my parameters
    potLoc = strippedSeq.find(forward)
    filterBool = (filterSeq(potLoc,forward) and filterSeq(-1,backward))
    primerLengthBool = inRange(len(forward),18,22) and inRange(len(backward),18,22)
    CGContentBool = inRange(pCG(forward),40,60) and inRange(pCG(backward),40,60)
    repeatsBool = fRepeats(forward) == -1 and fRepeats(backward) == -1
    tmBool = abs(getTmScore(forward)-getTmScore(backward)) <= 10
    print(filterBool,primerLengthBool,CGContentBool,repeatsBool,tmBool)

referencePrimerLocs = set()

all_alignments = []
index = 1

for (theirF, theirR) in referencePrimers:
    referencePrimerLocs.add((findLoc(theirF.lower(), False), findLoc(theirR.lower(), True)))

# To see if intervals overlap
def doesOverlap (s1,s2):
    (a,b) = s1
    (c,d) = s2
    return (c >= a and c <= b) or (d >= a and d <= b) or (a >= c and a <= d) or (b >= c and b <= d)

for (a,b) in candidatePrimerLocs:
    for (c,d) in referencePrimerLocs:
         all_alignments.append("Set " + str(index) + "\n")
         all_alignments.append("Forward: (" + str(a) + "," + str(c) + ") -> DELTA:" + str(a-c) +"\n")
         all_alignments.append("Backward: (" + str(b) + "," + str(d) + ") -> DELTA:" + str(b-d) +"\n")
         all_alignments.append(str(doesOverlap((a,b),(c,d))) + "\n")
         index += 1 

f = open("differenceOutputs-59.txt", "w") 
for a in all_alignments:
    f.write(a)
f.close() 