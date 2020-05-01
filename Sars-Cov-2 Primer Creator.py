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
for seq_record in SeqIO.parse("all-seqAlign-59.fasta", "fasta"):
    # print(seq_record.id)
    arr[i] = str(seq_record.seq).lower()
    strippedArr = arr[i].strip("-")
    # print(len(seq_record))
    i += 1

# For interested sequence
ignore = 0
intSeq = str(arr[ignore])
# stripped version
strippedSeq = intSeq.replace("-","")

''' Fix for issue with blanks in main sequence. 
        Need to create stripped version of strings, also map for them
            (And possibly clean the new interval)                    '''
reMap = dict()
counter = 0
for nuc in range(len(intSeq)):
    if intSeq[nuc] != "-":
        reMap[nuc] = counter
        counter += 1

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

# REQUIRES: 
#   - GC content parameter is satisfied
#   - Repeats parameter is satisfied
# ENSURES:
#   - Primer is reverse complement of string
#       starting at i, and being required length
#   - Primer length is 18-22 BP long
#   - Primer ends with G/C within last 5 bases
#       (but no more than 3 in that last 5)
def reversePrimer(start, end):
    five2Three = intSeq[start:(end+1)]
    compStr = genComplement(five2Three)
    # Takes complement (read 3' to 5' and flips to 5' to 3')
    return compStr[::-1]

# REQUIRES: 
#   - GC content parameter is satisfied
#   - Repeats parameter is satisfied
# ENSURES:
#   - Primer is complement of reverse complement of string
#       starting at i, and being required length
#   - Primer length is 18-22 BP long
#   - Primer ends with G/C within last 5 bases
#       (but no more than 3 in that last 5)
def forwardPrimer(start, end):
    five2Three = intSeq[start:(end+1)]
    # No need for complement since we get the complement of the complement,
    #   and resulting string is already 5' to 3'.
    return five2Three

# Returns true if the resulting character is considered "unique" in the
#   alignment (everything else is blank except character in interested seq)
def findBlank(i):
    for j in range(0, ignore):
        if arr[j][i] != "-":
            return False
    for j in range(ignore + 1,len(arr)):
        if arr[j][i] != "-":
            return False
    return intSeq[i] != '-'

# # A more robust version of the above (if blanks + mismatches is above threshold%)
# def findBlank(i):
#     # threshold
#     threshold = 100
#     # Refers to A,T,C,G,- respectively
#     counts = [0,0,0,0,0]
#     diccy = {'a':0, 't':1, 'c':2, 'g':3, '-':4}
#     reverseDiccy = {0:'a', 1:'t', 2:'c', 3:'g', 4:'-'}
#     for j in range(0, ignore):
#         if arr[j][i] in diccy:
#             counts[diccy[arr[j][i]]] += 1
#     for j in range(ignore + 1,len(arr)):
#         if arr[j][i] in diccy:
#             counts[diccy[arr[j][i]]] += 1
#     tot = sum(counts)
#     counts[diccy[intSeq[i]]] = 0
#     return sum(counts)/tot*100 >= threshold and intSeq[i] != '-'

# Provides CG content percentage in range [start,end] (note the [], not ())
def percentCG(start,end):
    ret = 0
    for i in range(start, end+1):
        if intSeq[i] == 'g' or intSeq[i] == 'c':
            ret += 1
    return ret/(end-start+1)*100

# Sliding window, if contents of window are 4-repeat we reveal it to be bad
def findRepeats(start,end):
    # If such a 4-long repeat exists in any form then we toss it out
    i = start
    s = set()
    s.add("aaaa")
    s.add("tttt")
    s.add("cccc")
    s.add("gggg")
    while (i+4 <= end):
        if intSeq[i:(i+4)] in s:
            return i
        i += 1
    return -1

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

# Forward Primer creation
def shiftStart(newStart, sShift, peak):
    while (newStart + sShift < peak):
        # Get out of repeats block, reset window
        repSLoc = findRepeats(newStart, newStart + sShift)
        if repSLoc != -1:
            newStart = repSLoc + 2
            sShift = 18
        # Get in the 40-60% CG range by increasing window size
        pCG = percentCG(newStart, newStart + sShift)
        # Use 45 to give wiggle room when increasing range
        #   Result can worst case lie just shy of 40% (like -2%)
        if pCG >= 45.0 and pCG <= 60.0:
            # Shift window to end with C/G (possibly)
            while (newStart + sShift < peak and sShift < 22):
                if intSeq[newStart+sShift] == 'c' or intSeq[newStart+sShift] == 'g':
                    giveUp = findRepeats(newStart, newStart + sShift)
                    if giveUp == -1:
                        return (newStart, sShift)
                    else: # The largest primer we can get
                        return (newStart, giveUp - newStart)
                sShift += 1
            # print(newStart, sShift, forwardPrimer(newStart, newStart+sShift))
        else:
            newStart += 1 # to potentially get into the percentage window
    return (-1, -1)

# Reverse Primer creation
def shiftEnd(newEnd, eShift, pit):
    while (newEnd - eShift > pit):
        # Get out of repeats block, reset window
        repELoc = findRepeats(newEnd - eShift, newEnd)
        if repELoc != -1:
            newEnd = repELoc - 1
            eShift = 18
        # Get in the 40-60% CG range by increasing window size
        pCG = percentCG(newEnd - eShift, newEnd)
        # Use 47 to give wiggle room when increasing range
        if pCG >= 47.0 and pCG <= 60.0:
            # Shift window to end with C/G (possibly)
            while (newEnd - eShift > pit and eShift < 22):
                if intSeq[newEnd - eShift] == 'c' or intSeq[newEnd - eShift] == 'g':
                    giveUp = findRepeats(newEnd - eShift, newEnd)
                    if giveUp == -1:
                        return (newEnd, eShift)
                    else: # The largest primer we can get
                        return (newEnd, newEnd - giveUp + 1)
                eShift += 1
        newEnd -= 1 # to potentially get into the percentage window
    return (-1, -1)

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
    maxScore = -1.0
    for i in range(len(strippedSeq)-len(seq)):
        # Calculate score
        tempScore = 0.0
        for j in range(len(seq)):
            if seq[j] == strippedSeq[i+j]:
                # Extreme bias in weighting for matches at 3' end, since last
                #   5 nucleotides at 3' end drastically affect how well primer
                #       anneals to strand
                if isRev == False and j > len(seq) - 5:
                    tempScore += 1 # Bias for forward primer
                if isRev == True and j < 5:
                    tempScore += 1 # Bias for reverse primer
                tempScore += 1.0
        if tempScore > maxScore:
            maxScore = tempScore
            maxIndex = i
    return maxIndex

# Splits a range that is larger than 500 into overlapping windows, and returns 
#   the range otherwise.
def splitRanges(begin, currNuc):
    dist = currNuc - begin
    acc = []
    if (dist > 500):
        acc = acc + splitRanges(begin, begin + dist//2)
        acc = acc + splitRanges(begin + dist//2, currNuc)
        acc = acc + splitRanges(begin + dist//4, begin + 3*dist//4)
    else:
        acc.append((begin, currNuc))
    return acc

# # Finds candidate ranges to begin looking
# def findCandidateRegions ():
#     threshold = 100
#     s = set()
#     begin = 0
#     tracking = False

#     for currNuc in range(len(intSeq)):
#         # Start looking for blanks (as in, a unique region in the aligned sequences)
#         if findBlank(currNuc):
#             if not tracking:
#                 # Mark "ranges" where I can begin to look. 
#                 tracking = True
#                 begin = currNuc
#         else:
#             if tracking:
#                 tracking = False
#                 dist = currNuc - begin
#                 # if the range is extremely large then cut in half
#                 #   also add overlapping middle range 
#                 acc = splitRanges(begin, currNuc)
#                 for (st, ed) in acc:
#                     # Keep the ones that exceed the threshold.
#                     if (ed - st) >= threshold:
#                         s.add((st, ed))
#     return s

# Finds candidate ranges to begin looking (Version 2)
def findCandidateRegions ():
    s = []
    begin = 0
    tracking = False
    ''' Parameters '''
    # Minimum size of amplified region
    threshold = 100
    # Minimum number of unique regions to be considered
    minUniques = 3

    # Create the set of unique regions
    for currNuc in range(len(intSeq)):
        # Start looking for blanks (as in, a unique region in the aligned sequences)
        if findBlank(currNuc):
            if not tracking:
                # Mark "ranges" where I can begin to look. 
                tracking = True
                begin = currNuc
        else:
            if tracking:
                tracking = False
                dist = currNuc - begin
                s.append((begin, currNuc))
    sorted(s, key=lambda x: x[0])
    retSet = set()
    i = 0
    sUpdated = []
    for (m,n) in s:
        sUpdated.append((reMap[m], reMap[n]))
    s = sUpdated
    while (i < len(s)):
        numUniques = i
        (a,b) = s[i]
        beg = a
        end = b
        while (b - beg < 500 and i < len(s)):
            end = b
            (a,b) = s[i]
            i += 1
        numUniques = i - numUniques
        if numUniques > minUniques:
            retSet.add((beg,end))
    return retSet

s = findCandidateRegions()
# need to strip any "-" from intseq and remap the stuff with it too 
intSeq = strippedSeq
print(s)

# Creates Primers
candidatePrimers = set()
candidatePrimerLocs = set()
for (start,end) in s:
    (newStart, sShift) = shiftStart(start, 18, end - 18)
    (newEnd, eShift) = shiftEnd(end, 18, newStart + sShift)
    # Primers of acceptable condition can't be created, or...
    if(newStart == -1 or newEnd == -1):
        continue 
    # Length of resulting amplified DNA is too short
    if(newEnd - newStart + 1 < 100):
        continue
    fPrimer = forwardPrimer(newStart, newStart + sShift)
    rPrimer = reversePrimer(newEnd - eShift, newEnd)
    # More filtering for TM difference
    while (abs(getTmScore(fPrimer) - getTmScore(rPrimer)) > 10):
        (newStart, sShift) = shiftStart(start, 18, end - 18)
        (newEnd, eShift) = shiftEnd(end, 18, start + 18)
        if(newStart == -1 or newEnd == -1):
            continue 
        fPrimer = forwardPrimer(newStart, newStart + sShift)
        rPrimer = reversePrimer(newEnd - eShift, newEnd)
    # Finally, check if the primer can be detected anywhere else
    if (filterSeq(newStart, fPrimer) and filterSeq(-1, rPrimer)):
        candidatePrimers.add((fPrimer, rPrimer))
        candidatePrimerLocs.add((newStart, newEnd-eShift))
print(candidatePrimers)

referencePrimers = set()

referencePrimers.add(("GGGTTGGGACTATCCTAAGTGTGA", "TAACACACAACACCATCATCA"))
referencePrimers.add(("AACACAAGCTTTCGGCAG", "GAAATTTGGATCTTTGTCATCC"))
referencePrimers.add(("TGCGGCCAATGTTTGTAATCAGCCAAGGAAATTTTGGGGAC", "CGCATTGGCATGGAAGTCACTTTGATGGCACCTGTGTAG"))
referencePrimers.add(("TTCCTTGTCTGATTAGTTC", "ACCTTCGGGAACGTGGTT"))


''' Tests for the reference primers
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
    filterBool = (filterSeq(potLoc,forward) and filterSeq(-1,forward))
    primerLengthBool = inRange(len(forward),18,22) and inRange(len(backward),18,22)
    print(pCG(forward), pCG(backward))
    CGContentBool = inRange(pCG(forward),40,60) and inRange(pCG(backward),40,60)
    repeatsBool = fRepeats(forward) == -1 and fRepeats(backward) == -1
    tmBool = abs(getTmScore(forward)-getTmScore(backward)) <= 10
    print(filterBool,primerLengthBool,CGContentBool,repeatsBool,tmBool)
 '''


referencePrimerLocs = set()

all_alignments = []
index = 1

for (theirF, theirR) in referencePrimers:
    referencePrimerLocs.add((findLoc(theirF.lower(), False), findLoc(theirR.lower(), True)))

for (a,b) in candidatePrimerLocs:
    for (c,d) in referencePrimerLocs:
         all_alignments.append("Set " + str(index) + "\n")
         all_alignments.append("Forward: (" + str(a) + "," + str(c) + ") -> DELTA:" + str(a-c) +"\n")
         all_alignments.append("Backward: (" + str(b) + "," + str(d) + ") -> DELTA:" + str(b-d) +"\n")
         index += 1 

f = open("alignedOutputs.txt", "w") 
for a in all_alignments:
    f.write(a)
f.close() 