#!/usr/bin/env python3
""""
@Name: find_repeats.py
@Author: mmcguffi

Creates a least-common-prefix array from a string and identifies repeats
@Todo: Improve docstring
"""
import argparse
from collections import defaultdict, ChainMap, namedtuple
from operator import itemgetter
from itertools import chain

# -- Public Facing Methods and Classes --

class RepeatAnalyzer():
    """Identifies repeated sequences and provides functions to retrieve them"""
    def __init__(self, sequence, minlength=0, filtered=True):
        # Set up variables
        self.sequence: str = str(sequence)
        self._repeat_db = {}
        self._process_palandromes(minlength, filtered)
        self.repeated_sequences: list = self._repeat_db.keys()

    def _process_palandromes(self, minlength, filtered):
        # Generate potential repeats and palandrome/hairpins
        table = self.sequence.maketrans("ATCGatcg", "TAGCtagc")
        reverse_seq = self.sequence.translate(table)
        concat_seq = self.sequence + reverse_seq
        sequence_length = len(self.sequence)
        if filtered is False:
            candidate_seqs = get_repeats(concat_seq, minlength)
        else:
            candidate_seqs = get_repeats_filtered(concat_seq, minlength)

        #Colapse length
        candidate_seqs = [candidate_seqs[key] for key in candidate_seqs]
        candidate_seqs = dict(ChainMap(*candidate_seqs))

        # Separate out palandromes from hairpins

        for candidate, locations in candidate_seqs.items():
            hairpin_positions = []
            palandrome_positions = []
            candidate_length = len(candidate)
            repeat_positions = [location for location in locations if location < sequence_length-candidate_length]

            potential_hairpins = [location for location in locations if location not in repeat_positions]
            potential_hairpins = [location-sequence_length for location in potential_hairpins
                                  if location-sequence_length >= 0]

            # Identifies palandromes and picks them out. Also does so if they are palandromes too.
            for position in repeat_positions:
                for potential_palandrome in potential_hairpins:
                    if position + candidate_length == potential_palandrome:
                        palandrome_positions.append(potential_palandrome)
                    if potential_palandrome + candidate_length < position or \
                                   position + candidate_length < potential_palandrome:
                        hairpin_positions.append(potential_palandrome)

            self._repeat_db[candidate] = {'repeats': set(repeat_positions),
                                          'palandromes': set(palandrome_positions),
                                          'hairpins': set(hairpin_positions)}

    def get_repeats(self, checked_seq: str) -> list:
        """Recieves string and returns where string is repeated in a sequence"""
        if checked_seq in self._repeat_db:
            return self._repeat_db[checked_seq]['repeats']
        return []

    def get_palandromes(self, checked_seq: str) -> list:
        """Recieves string and returns where palandromes begin"""
        if checked_seq in self._repeat_db:
            palandrome_middle = self._repeat_db[checked_seq]['palandromes']
            palandrome_start = [position-len(checked_seq) for position in palandrome_middle]
            return palandrome_start
        return []

    def get_hairpins(self, checked_seq: str) -> list:
        """Recieves string and returns the start of where it may form a hairpin"""
        if checked_seq in self._repeat_db:
            positions = self.get_repeats(checked_seq)
            rc_positions = self._repeat_db[checked_seq]['hairpins']
            Hairpin = namedtuple('Hairpin', ['stemf', 'stemr', 'loop'])
            length_checked_seq = len(checked_seq)
            hairpins = []
            already_registered = []
            for forward_pos in positions:
                for reverse_pos in rc_positions:
                    if forward_pos > reverse_pos and forward_pos + length_checked_seq - reverse_pos >= 6:
                        loop_length = forward_pos - length_checked_seq - reverse_pos
                        if [reverse_pos+length_checked_seq+1, reverse_pos+length_checked_seq+loop_length] in already_registered:
                            continue
                        hairpins.append(Hairpin([forward_pos, forward_pos+length_checked_seq],
                                                [reverse_pos, reverse_pos+length_checked_seq],
                                                [reverse_pos+length_checked_seq+1, reverse_pos+length_checked_seq+loop_length]))
                        already_registered.append([reverse_pos+length_checked_seq+1, reverse_pos+length_checked_seq+loop_length])
                    elif reverse_pos > forward_pos and reverse_pos + length_checked_seq - forward_pos >= 6:
                        loop_length = reverse_pos - length_checked_seq - forward_pos
                        if [forward_pos+length_checked_seq+1, forward_pos+length_checked_seq+loop_length] in already_registered:
                            continue
                        hairpins.append(Hairpin([forward_pos, forward_pos+length_checked_seq],
                                                [reverse_pos, reverse_pos+length_checked_seq],
                                                [forward_pos+length_checked_seq+1, forward_pos+length_checked_seq+loop_length]))
                        already_registered.append([forward_pos+length_checked_seq+1, forward_pos+length_checked_seq+loop_length])
            return hairpins
        return []

    def fetch_repeats_of_length(self,checked_length: int,
                                find_higher: bool = False) -> list:
        """Recieves length and returns repeats of that length"""
        if find_higher is False:
            valid_keys = [key for key in self._repeat_db if len(key) == checked_length]
        else:
            valid_keys = [key for key in self._repeat_db if len(key) >= checked_length]
        return valid_keys

    def fetch_repeats_with_freq(self, checked_number: int,
                                find_higher: bool = True,
                                include_rc: bool = False) -> list:
        """Recieves frequency and returns repeats of (at least) that frequency"""
        valid_keys = []
        if include_rc is True and find_higher is False:
            valid_keys = [repeat for repeat, positions in self._repeat_db.items()
                          if len(chain(*positions.values()) == checked_number)]
        elif include_rc is True and find_higher is True:
            valid_keys = [repeat for repeat, positions in self._repeat_db.items()
                          if len(chain(*positions.values()) >= checked_number)]
        elif include_rc is False and find_higher is False:
            valid_keys = [repeat for repeat, positions in self._repeat_db.items()
                          if len(chain(positions[0]) == checked_number)]
        elif include_rc is False and find_higher is True:
            valid_keys = [repeat for repeat, positions in self._repeat_db.items()
                          if len(chain(positions[0]) >= checked_number)]
        else:
            raise ValueError("Invalid Arguments")
        if not valid_keys:
            return []
        return valid_keys

def get_repeats(in_str, repeat_length=6):
    """""" # @TODO: Docstring
    lcparray = create_lcparray(in_str, repeat_length)
    repeat_dict = defaultdict(lambda: {})
    for i, _ in enumerate(lcparray):
        match = in_str[lcparray[i][0]:lcparray[i][0] + lcparray[i][1]]
        if match in repeat_dict[len(match)].keys():
            repeat_dict[len(match)][match].append(lcparray[i][0])
            repeat_dict[len(match)][match].append(lcparray[i][2])
        else:
            repeat_dict[len(match)][match] = [lcparray[i][0], lcparray[i][2]]
        repeat_dict[len(match)][match] = list(set(repeat_dict[len(match)][match]))
    return repeat_dict

def get_repeats_filtered(in_seq, repeat_length=6, return_dir_rep=False):
    """""" # @TODO: Docstring
    # essesntially creates a 1D sequence space denoted where repeats are, and uses set theory to
    ##check if a smaller repeat should be included
    repeat_space = [[] for i in range(len(in_seq))]
    unfiltered_repeats = get_repeats(in_seq, 1)
    direct_repeats = get_direct_repeats(unfiltered_repeats)
    unfiltered_repeats = {k: v for k, v in unfiltered_repeats.items() if k >= repeat_length}

    for len_rep, value in direct_repeats.items():
        for i, _ in enumerate(value):
            start = direct_repeats[len_rep][i][0]
            end = direct_repeats[len_rep][i][1]
            for j in range(end - start):
                repeat_space[start + j].append("dr" + str(len_rep) + "." + str(i))

    nested_dict = lambda: defaultdict(nested_dict)  # this allows for the generation of non-declared nested dicts
    filteredrepeat_dict = nested_dict()

    key_list = sorted(list(unfiltered_repeats.keys()), reverse=True)  # sorts keys highest to lowest
    counter = 0
    for key in key_list:
        for seqkey in unfiltered_repeats[key]:
            counter += 1
            set_space = []
            if not check_overlap(unfiltered_repeats[key][seqkey], key):  # filters out "non-anchored" overlaps
                for pos in unfiltered_repeats[key][seqkey]:
                    set_space.append(set_intersection(repeat_space[pos:pos + key]))
                condensedset_space = set_intersection(set_space)
                if check_clean_dr_space(condensedset_space):
                    if len(condensedset_space) == 0:
                        for pos in unfiltered_repeats[key][seqkey]:
                            filteredrepeat_dict[key][seqkey] = unfiltered_repeats[key][seqkey]
                            for i in range(key):
                                repeat_space[pos + i].append(counter)
    if return_dir_rep is False:
        return default_to_regular(filteredrepeat_dict)
    if return_dir_rep is True:
        return default_to_regular(filteredrepeat_dict), direct_repeats  # , repeat_space
    raise ValueError("return_dir_rep must be boolian")


# basic repeat methods
def create_lcparray(in_str, repeat_length):
    """""" # @TODO: Docstring
    suffixarray = []
    for i in range(1, len(in_str) + 1):
        suffixarray.append((len(in_str) - i, in_str[-i:len(in_str)]))
    suffixarray.sort(key=itemgetter(1))

    lcparray = []
    for i in range(len(suffixarray) - 1):
        j = 0
        try:
            while suffixarray[i][1][j] == suffixarray[i + 1][1][j]:
                j += 1
                if j >= repeat_length:
                    lcparray.append((suffixarray[i][0], j, suffixarray[i + 1][0]))
        except IndexError:
            pass

    lcparray.sort(key=itemgetter(1), reverse=True)

    return lcparray

# advanced repeat methods

def default_to_regular(defdict: defaultdict) -> dict:
    """Recursively converts input defaultdict to a returned dict"""
    if isinstance(defdict, defaultdict):
        defdict = {k: default_to_regular(v) for k, v in defdict.items()}
    return defdict


def set_intersection(list_of_sublists: list) -> list:
    """Returns list of elements common to all sublists"""
    result = set(list_of_sublists[0])
    for sublist in list_of_sublists[1:]:
        result.intersection_update(sublist)
    return list(result)


def check_clean_dr_space(in_condensed_set_space):
    """""" # @TODO: Multiline docstring -> This is technical
    # checks if seq comparison is entirely within a direct repeat space
    for i in in_condensed_set_space:
        if 'dr' in str(i):
            return False
    return True


def check_overlap(in_pos_list, in_len):
    """""" # @TODO: Multiline docstring -> This is technical
    # if False, it means there is either no overlap, or it is "anchored" by a third
    # if True, it means they are overlapping
    in_pos_list.sort()
    shifted_in_pos_list = [x + in_len for x in in_pos_list]
    for i in range(len(in_pos_list) - 1):
        if in_pos_list[i + 1] >= shifted_in_pos_list[i]:
            return False
    return True


def get_direct_repeats(in_dict):
    """""" # @TODO: Docstring
    directrepeat_dict = defaultdict(lambda: [])
    for seqLen in in_dict:
        for key in in_dict[seqLen]:
            if (key + key)[1:-1].find(key) != -1: ##checks to see if it can be broken down into smaller repeating substrings
                continue  ##if it is a unique seq, false -- if it repeats seq, true
            key_len = len(key)
            in_dict[seqLen][key].sort()
            j=0
            for i in range(len(in_dict[seqLen][key])-1):
                if j != 0: #this "remembers" where we were in the loop, so we dont re-count the same repeats
                    j-=1 
                else:
                    if in_dict[seqLen][key][i+j]+len(key) == in_dict[seqLen][key][i+j+1]:
                        try:
                            while in_dict[seqLen][key][i+j]+len(key) == in_dict[seqLen][key][i+j+1]:
                                j+=1
                        except IndexError:
                            pass

                        #foundSeq = get_feature_seq(ecoli.features[CDS])[inDict[seqLen][key][i]:inDict[seqLen][key][i]+((j+1)*len(key))]

                        if (key_len == 1 and j + 1 >= 5) or\
                            (key_len == 2 and j + 1 >= 4) or\
                            (key_len < 15 and j + 1 >= 3) or\
                            (key_len >= 15 and j + 1 >= 2):
                            # print("single","\t",repeat_dict[seq_len][key][i],"\t",foundSeq)
                            value = in_dict[seqLen][key]
                            directrepeat_dict[key_len].append((value[i], value[i] + (j + 1) * key_len))
    filtered_direct_repeat_dict = {}
    # homopolyList=[]
    for key, value in directrepeat_dict.items():
        if key == 1:  # already "actively filtered" above
            filtered_direct_repeat_dict[key] = value
        else:
            filtered_direct_repeat_dict[key] = []
            value.sort(key=itemgetter(0))
            i = 0
            while i < len(value):
                j = 1
                try:
                    while value[i][0] + key - 1 >= value[i + j][0]:
                        j += 1
                except IndexError:
                    pass
                filtered_direct_repeat_dict[key].append(value[i])
                i += j

    return filtered_direct_repeat_dict



if __name__ == "__main__":
    # Test case
    M_CHERRY = "ATGGTATCTAAGGGTGAAGAAGACAACATGGCTATCATCAAGGAATTCATGCGCTTCAAAGTTCACATGGAAGGCTCCGTTAACGGCCACGAATTCGAAATCG\
                AAGGCGAAGGCGAAGGTCGTCCGTACGAAGGTACTCAGACCGCTAAACTGAAAGTAACTAAAGGCGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCA\
                GTTCATGTACGGTTCCAAAGCTTACGTAAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCTTTCCCGGAAGGCTTCAAATGGGAACGTGTAATGAAC\
                TTCGAAGACGGCGGTGTTGTTACCGTTACTCAGGACTCTTCCCTGCAGGACGGCGAATTCATCTACAAGGTAAAGCTGCGCGGTACCAACTTCCCGTCCGACG\
                GCCCGGTAATGCAGAAGAAAACTATGGGCTGGGAAGCTTCCTCTGAACGTATGTACCCGGAAGACGGCGCTCTGAAGGGTGAAATCAAACAGCGCCTGAAGCT\
                GAAGGACGGCGGCCACTACGACGCTGAAGTTAAAACCACTTACAAAGCTAAAAAACCGGTTCAGCTGCCGGGTGCTTACAACGTAAACATCAAGCTGGACATC\
                ACTTCCCACAACGAAGACTACACCATCGTAGAACAGTACGAACGTGCTGAAGGTCGCCACTCTACTGGCGGTATGGACGAACTGTACAAG"
    M_CHERRY = RepeatAnalyzer(M_CHERRY)
    print(M_CHERRY.repeated_sequences)


    parser = argparse.ArgumentParser(description="find repeats in strs")
    parser.add_argument("-f", "--file_path", action="store", help="path to txt file containing only DNA bases")
    parser.add_argument("-l", "--length", action="store", help="length of repeats to look for", default=6)
    parser.add_argument("-d", "--direct_repeats", action="store_true", help="also find direct repeats", default=False)
 
    args = parser.parse_args()
    """
    with open(args.file_path, 'r', encoding='utf8') as f:
        seq = f.read()
    seq = seq.replace('\n', '')
    seq = seq.replace('\r', '')
    seq = seq.replace(' ', '')
    seq = seq.upper()

    hits = get_repeats_filtered(seq, int(args.length), args.direct_repeats)
    print(hits)"""
# EOF
