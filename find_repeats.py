#!/usr/bin/env python3
import argparse
from collections import defaultdict
from operator import itemgetter


#basic repeat methods
def create_lcparray(inStr, repeat_length):
    suffixarray = []
    for i in range(1,len(inStr)+1):
        suffixarray.append((len(inStr)-i,inStr[-i:len(inStr)]))
    suffixarray.sort(key=itemgetter(1))

    lcparray=[]
    for i in range(len(suffixarray)-1):
        j=0
        try:
            while suffixarray[i][1][j]==suffixarray[i+1][1][j]:
                j+=1
                if j>=repeat_length:
                    lcparray.append((suffixarray[i][0],j,suffixarray[i+1][0]))
        except IndexError:
            pass

    lcparray.sort(key=itemgetter(1), reverse=True)

    return lcparray
def get_repeats(inStr, repeat_length = 6):
    lcparray = create_lcparray(inStr, repeat_length)
    repeatDict = {}
    for i in range(len(lcparray)):
        match = inStr[lcparray[i][0]:lcparray[i][0]+lcparray[i][1]]
        try:
            repeatDict[len(match)]
        except KeyError:
            repeatDict[len(match)]={}
        try:
            repeatDict[len(match)][match].append(lcparray[i][0])
            repeatDict[len(match)][match].append(lcparray[i][2])
        except:
            repeatDict[len(match)][match] =[lcparray[i][0],lcparray[i][2]]

        repeatDict[len(match)][match] = list(set(repeatDict[len(match)][match]))
        
    return repeatDict


#advanced repeat methods
def get_repeats_filtered(inSeq,repeat_length = 6,returnDirRep=False):
    #essesntially creates a 1D sequence space denoted where repeats are, and uses set theory to
    ##check if a smaller repeat should be included
    repeatSpace= [[] for i in range(len(inSeq))]
    unfilteredRepeats = get_repeats(inSeq,1)
    directRepeats=get_direct_repeats(unfilteredRepeats)
    unfilteredRepeats={k: v for k, v in unfilteredRepeats.items() if k>=repeat_length}
    
    for lenRep in directRepeats:
        for i in range(len(directRepeats[lenRep])):
            start=directRepeats[lenRep][i][0]
            end=directRepeats[lenRep][i][1]
            for j in range(end-start):
                repeatSpace[start+j].append("dr"+str(lenRep)+"."+str(i))

    nested_dict = lambda: defaultdict(nested_dict) #this allows for the generation of non-declared nested dicts
    filteredRepeatDict = nested_dict()
    
    keyList=sorted(list(unfilteredRepeats.keys()),reverse=True) #sorts keys highest to lowest
    c=0
    for key in keyList:
        for seqkey in unfilteredRepeats[key]:
            c += 1
            setSpace=[]
            if not check_overlap(unfilteredRepeats[key][seqkey],key): #filters out "non-anchored" overlaps
                for pos in unfilteredRepeats[key][seqkey]:
                    setSpace.append(set_intersection(repeatSpace[pos:pos+key]))
                condensedSetSpace=set_intersection(setSpace)
                if check_clean_dr_space(condensedSetSpace):
                    if len(condensedSetSpace) == 0:
                         for pos in unfilteredRepeats[key][seqkey]:
                            filteredRepeatDict[key][seqkey]=unfilteredRepeats[key][seqkey]
                            for i in range(key):
                                repeatSpace[pos+i].append(c)
    if returnDirRep == False:
        return default_to_regular(filteredRepeatDict)
    elif returnDirRep == True:
        return default_to_regular(filteredRepeatDict), directRepeats#, repeatSpace
    else:
        print("error")

def default_to_regular(d): #default dict to reg dict, recusive for arb dict depth
    if isinstance(d, defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d
def set_intersection(listOfSublists):
    #takes a list of sublists and returns a list of all common elements in sublist
    result = set(listOfSublists[0])
    for s in listOfSublists[1:]:
        result.intersection_update(s)
    return list(result)
def check_clean_dr_space(inCondensedSetSpace):
    #checks if seq comparison is entirely within a direct repeat space
    for i in inCondensedSetSpace:
        if 'dr' in str(i):
            return False
    return True
def check_overlap(inPosList,inLen):
    #if False, it means there is either no overlap, or it is "anchored" by a third
    #if True, it means they are overlapping
    inPosList.sort()
    shiftedInPosList=[x+inLen for x in inPosList]
    for i in range(len(inPosList)-1):
        if inPosList[i+1] >= shiftedInPosList[i]:
            return False
    return True

def get_direct_repeats(inDict):
    directrepeatDict={}

    for seqLen in inDict:
        for key in inDict[seqLen]:
            if (key + key)[1:-1].find(key) != -1: ##checks to see if it can be broken down into smaller repeating substrings
                pass  ##if it is a unique seq, false -- if it repeats seq, true
            else:
                inDict[seqLen][key].sort()
                j=0
                for i in range(len(inDict[seqLen][key])-1):
                    if j != 0: #this "remembers" where we were in the loop, so we dont re-count the same repeats
                        j-=1 
                        pass
                    else:
                        if inDict[seqLen][key][i+j]+len(key) == inDict[seqLen][key][i+j+1]:
                            try:
                                while inDict[seqLen][key][i+j]+len(key) == inDict[seqLen][key][i+j+1]:
                                    j+=1
                            except IndexError:
                                pass

                            #foundSeq = get_feature_seq(ecoli.features[CDS])[inDict[seqLen][key][i]:inDict[seqLen][key][i]+((j+1)*len(key))]

                            if len(key)==1:
                                if j+1 >= 5:
                                    #print("single","\t",repeatDict[seqLen][key][i],"\t",foundSeq)
                                    try:
                                        directrepeatDict[len(key)].append((inDict[seqLen][key][i],inDict[seqLen][key][i]+(j+1)*len(key)))
                                    except KeyError:
                                        directrepeatDict[len(key)]=[(inDict[seqLen][key][i],inDict[seqLen][key][i]+(j+1)*len(key))]
                            elif len(key)==2:
                                if j+1 >= 4:
                                    #print("direp","\t",repeatDict[seqLen][key][i],"\t",foundSeq)
                                    try:
                                        directrepeatDict[len(key)].append((inDict[seqLen][key][i],inDict[seqLen][key][i]+(j+1)*len(key)))
                                    except KeyError:
                                        directrepeatDict[len(key)]=[(inDict[seqLen][key][i],inDict[seqLen][key][i]+(j+1)*len(key))]
                            elif len(key)<15:
                                if j+1 >= 3:
                                    #print(str(seqLen)+"mer","\t",repeatDict[seqLen][key][i],"\t",foundSeq)
                                    try:
                                        directrepeatDict[len(key)].append((inDict[seqLen][key][i],inDict[seqLen][key][i]+(j+1)*len(key)))
                                    except KeyError:
                                        directrepeatDict[len(key)]=[(inDict[seqLen][key][i],inDict[seqLen][key][i]+(j+1)*len(key))]
                            elif len(key)>=15:
                                if j+1 >= 2:
                                    #print(str(seqLen)+"mer","\t",repeatDict[seqLen][key][i],"\t",foundSeq)
                                    try:
                                        directrepeatDict[len(key)].append((inDict[seqLen][key][i],inDict[seqLen][key][i]+(j+1)*len(key)))
                                    except KeyError:
                                        directrepeatDict[len(key)]=[(inDict[seqLen][key][i],inDict[seqLen][key][i]+(j+1)*len(key))]
    filteredDirectrepeatDict={}
    #homopolyList=[]
    for key in directrepeatDict:
        if key == 1: #already "actively filtered" above
            filteredDirectrepeatDict[key]=directrepeatDict[1]
        else:
            filteredDirectrepeatDict[key]=[]
            directrepeatDict[key].sort(key=itemgetter(0))
            i=0
            while i < len(directrepeatDict[key]):
                j=1
                try:
                    while directrepeatDict[key][i][0]+key-1>=directrepeatDict[key][i+j][0]:
                        j+=1
                except IndexError:
                    pass
                filteredDirectrepeatDict[key].append(directrepeatDict[key][i])
                i+=j
    
    return filteredDirectrepeatDict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="find repeats in strs")
    parser.add_argument("-f", "--file_path", action="store", help="path to txt file containing only DNA bases")
    parser.add_argument("-l", "--length", action="store", help="length of repeats to look for" , default=6)
    parser.add_argument("-d", "--direct_repeats", action="store_true", help="also find direct repeats", default=False)


    args = parser.parse_args()
    
    with open(args.file_path, 'r') as f:
        seq = f.read()
    seq = seq.replace('\n', '')
    seq = seq.replace('\r', '')
    seq = seq.replace(' ', '')
    seq = seq.upper()

    hits = get_repeats_filtered(seq, int(args.length), args.direct_repeats)
    print(hits)
