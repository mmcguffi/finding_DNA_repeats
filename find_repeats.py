#!/usr/bin/env python3
""""
@Name: find_repeats.py
@Author: mmcguffi

Creates a least-common-prefix array from a string and identifies repeats
@Todo: Improve docstring
"""

import argparse
from collections import defaultdict
from operator import itemgetter


# --- User-facing functions ---
def get_repeats(instr: str, repeat_length: int = 6) -> dict:
    """"""  # @Todo: Add docstring
    lcparray = create_lcparray(instr, repeat_length)
    repeatDict = {}
    for i in range(len(lcparray)):
        match = instr[lcparray[i][0]:lcparray[i][0] + lcparray[i][1]]
        if len(match) not in repeatDict.keys():
            repeatDict[len(match)] = {}
        if match in repeatDict[len(match)].keys():
            repeatDict[len(match)][match].append(lcparray[i][0])
            repeatDict[len(match)][match].append(lcparray[i][2])
        else:
            repeatDict[len(match)][match] = [lcparray[i][0], lcparray[i][2]]

        repeatDict[len(match)][match] = list(set(repeatDict[len(match)][match]))
    return repeatDict


def get_repeats_filtered(instr: str, repeat_length: int = 6, return_dir_rep: bool = False):
    """Essesntially creates a 1D sequence space denoted where repeats are, and uses set theory to
    check if a smaller repeat should be included"""  # @Todo: Improve
    repeatSpace = [[] for i in range(len(instr))]
    unfiltered_repeats = get_repeats(instr, 1)
    direct_repeats = get_direct_repeats(unfiltered_repeats)
    unfiltered_repeats = {k: v for k, v in unfiltered_repeats.items() if k >= repeat_length}

    for lenRep in direct_repeats:
        for i in range(len(direct_repeats[lenRep])):
            start = direct_repeats[lenRep][i][0]
            end = direct_repeats[lenRep][i][1]
            for j in range(end - start):
                repeatSpace[start + j].append(["dr" + str(lenRep) + "." + str(i)])

    nested_dict = lambda: defaultdict(nested_dict)  # this allows for the generation of non-declared nested dicts
    filtered_repeat_dict = nested_dict()

    key_list = sorted(list(unfiltered_repeats.keys()), reverse=True)  # sorts keys highest to lowest
    c = 0
    for key in key_list:
        for seqkey in unfiltered_repeats[key]:
            c += 1
            set_space = []
            if not check_overlap(unfiltered_repeats[key][seqkey], key):  # filters out "non-anchored" overlaps
                for pos in unfiltered_repeats[key][seqkey]:
                    set_space.append(set_intersection(repeatSpace[pos:pos + key]))
                condensed_set_space = set_intersection(set_space)
                if check_clean_dr_space(condensed_set_space):
                    if len(condensed_set_space) == 0:
                        for pos in unfiltered_repeats[key][seqkey]:
                            filtered_repeat_dict[key][seqkey] = unfiltered_repeats[key][seqkey]
                            for i in range(key):
                                repeatSpace[pos + i].append([c])
    if return_dir_rep is False:
        return default_to_regular(filtered_repeat_dict), None
    elif return_dir_rep is True:
        return default_to_regular(filtered_repeat_dict), direct_repeats  # , repeatSpace


# --- Utility functions ---

#basic repeat methods
def create_lcparray(instr: str, repeat_length: int):
    """"""  # @Todo: Add docstring
    suffixarray = []
    for i in range(1,len(instr)+1):
        suffixarray.append((len(instr)-i,instr[-i:len(instr)]))
    suffixarray.sort(key=itemgetter(1))

    lcparray = []
    for i in range(len(suffixarray)-1):
        j = 0
        try:
            while suffixarray[i][1][j] == suffixarray[i+1][1][j]:
                j += 1
                if j >= repeat_length:
                    lcparray.append((suffixarray[i][0], j, suffixarray[i+1][0]))
        except IndexError:
            pass

    lcparray.sort(key=itemgetter(1), reverse=True)

    return lcparray


# advanced repeat methods
def default_to_regular(d: dict) -> dict:  # default dict to reg dict, recusive for arb dict depth
    """"""  # @Todo: Add docstring
    if isinstance(d, defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d


def set_intersection(sublists: list) -> list:
    """"""  # @Todo: Add docstring
    # takes a list of sublists and returns a list of all common elements in sublist
    result = set(sublists[0])
    for s in sublists[1:]:
        result.intersection_update(s)
    return list(result)


def check_clean_dr_space(in_condensed_set_space) -> bool:
    """"""  # @Todo: Add docstring
    #checks if seq comparison is entirely within a direct repeat space
    for i in in_condensed_set_space:
        if 'dr' in str(i):
            return False
    return True


def check_overlap(in_pos_list, in_len) -> bool:
    """"""  # @Todo: Add docstring
    # if False, it means there is either no overlap, or it is "anchored" by a third
    # if True, it means they are overlapping
    in_pos_list.sort()
    shifted_in_pos_list = [x+in_len for x in in_pos_list]
    for i in range(len(in_pos_list)-1):
        if in_pos_list[i+1] >= shifted_in_pos_list[i]:
            return False
    return True


def get_direct_repeats(in_dict: dict) -> dict:
    """"""  # @Todo: Add docstring
    directrepeats = {}
    for seqLen in in_dict:
        for key in in_dict[seqLen]:
            if (key + key)[1:-1].find(key) != -1: # checks to see if it can be broken down into smaller repeating substrings
                continue  # if it is a unique seq, false -- if it repeats seq, true
            else:
                in_dict[seqLen][key].sort()
                j = 0
                for i in range(len(in_dict[seqLen][key])-1):
                    if j != 0:  # this "remembers" where we were in the loop, so we dont re-count the same repeats
                        j -= 1
                    else:
                        if in_dict[seqLen][key][i+j]+len(key) == in_dict[seqLen][key][i+j+1]:
                            try:
                                while in_dict[seqLen][key][i+j]+len(key) == in_dict[seqLen][key][i+j+1]:
                                    j += 1
                            except IndexError:
                                pass

                            #foundSeq = get_feature_seq(ecoli.features[CDS])[in_dict[seqLen][key][i]:in_dict[seqLen][key][i]+((j+1)*len(key))]

                            if len(key) == 1:
                                if j+1 >= 5:
                                    #print("single","\t",repeatDict[seqLen][key][i],"\t",foundSeq)
                                    try:
                                        directrepeats[len(key)].append((in_dict[seqLen][key][i],
                                                                        in_dict[seqLen][key][i] + (j + 1) * len(key)))
                                    except KeyError:
                                        directrepeats[len(key)] = [(in_dict[seqLen][key][i],
                                                                    in_dict[seqLen][key][i] + (j + 1) * len(key))]
                            elif len(key) == 2:
                                if j+1 >= 4:
                                    #print("direp","\t",repeatDict[seqLen][key][i],"\t",foundSeq)
                                    if len(key) in directrepeats.keys():
                                        directrepeats[len(key)].append((in_dict[seqLen][key][i],
                                                                        in_dict[seqLen][key][i] + (j + 1) * len(key)))
                                    else:
                                        directrepeats[len(key)] = [(in_dict[seqLen][key][i],
                                                                    in_dict[seqLen][key][i] + (j + 1) * len(key))]
                            elif len(key) < 15:
                                if j+1 >= 3:
                                    #print(str(seqLen)+"mer","\t",repeatDict[seqLen][key][i],"\t",foundSeq)
                                    if len(key) in directrepeats.keys():
                                        directrepeats[len(key)].append((in_dict[seqLen][key][i],
                                                                        in_dict[seqLen][key][i] + (j + 1) * len(key)))
                                    else:
                                        directrepeats[len(key)] = [(in_dict[seqLen][key][i],
                                                                    in_dict[seqLen][key][i] + (j + 1) * len(key))]
                            elif len(key) >= 15:
                                if j+1 >= 2:
                                    #print(str(seqLen)+"mer","\t",repeatDict[seqLen][key][i],"\t",foundSeq)
                                    if len(key) in directrepeats.keys():
                                        directrepeats[len(key)].append((in_dict[seqLen][key][i],
                                                                        in_dict[seqLen][key][i] + (j + 1) * len(key)))
                                    else:
                                        directrepeats[len(key)] = [(in_dict[seqLen][key][i],
                                                                    in_dict[seqLen][key][i] + (j + 1) * len(key))]
    filtereddirectrepeats={}
    #homopolyList=[]
    for key in directrepeats:
        if key == 1:  # already "actively filtered" above
            filtereddirectrepeats[key] = directrepeats[1]
            continue
        else:
            filtereddirectrepeats[key]=[]
            directrepeats[key].sort(key=itemgetter(0))
            i = 0
            while i < len(directrepeats[key]):
                j = 1
                try:
                    while directrepeats[key][i][0]+key-1 >= directrepeats[key][i+j][0]:
                        j += 1
                except IndexError:
                    pass
                filtereddirectrepeats[key].append(directrepeats[key][i])
                i += j
    
    return filtereddirectrepeats


# ---- Run as script/Argparse ---
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
