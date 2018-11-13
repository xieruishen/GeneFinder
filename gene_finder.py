# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Sherrie Shen

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    """
    # TODO: implement this
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    pass

def get_reverse(dna):
    """ Computes the reverse sequence of DNA for the specified DNA sequences

        dna: a DNA sequence represented as a string
        returns: the reverse DNA sequence represented as a string
    >>> get_reverse("ATACCTGACGGT")
    'TGGCAGTCCATA'
    """
    new_dna = dna[::-1]
    return new_dna
    pass

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    new_dna = get_reverse(dna)
    reverse_complement = ''
    for letter in new_dna:
        complement = get_complement(letter)
        reverse_complement += complement

    return reverse_complement
    pass

def find_stop_codon(dna):
    """Takes a DNA sequence and search for a stop codon 'TAG' or 'TAA' or 'TGA'.
    If there is a stop codon, returns index of first letter of stop codon.
    If there is no stop codon, return index of the last letter in the sequence.
    >>> find_stop_codon("ATTCTAGCC")
    4
    >>> find_stop_codon("ATTCTACTAAGCC")
    7
    >>> find_stop_codon("CTACTGGATGACC")
    8
    >>> find_stop_codon("CCTCGGCTACG")
    11
    """
    index = 0
    while index < len(dna):
        if dna[index:index+3] == 'TAG':
            return index
        elif dna[index:index+3] == 'TAA':
            return index
        elif dna[index:index+3] == 'TGA':
            return index
        index = index + 1
    return len(dna)
    pass


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: implement this
    index = dna.find('ATG')
    #print(index)
    newindex = index+3
    stop_index = find_stop_codon(dna[newindex:]) + newindex
    #print(dna[newindex:])
    #print(stop_index)
    ORF = dna[index:stop_index]
    return ORF

    pass

def find_start_codon_3(dna):
    """ Return the index of the start codon in the sequence. If there is no start codon, return None.

      dna = a DNA sequences
      returns: The index of the start codon if there is one. None if there is none.
    >>> find_start_codon_3('ATAATGAGG')
    3
    """

    index = 0
    while index < len(dna):
        if dna[index:index+3] == 'ATG':
            return index
        else:
            index += 3
    return None


def find_stop_codon_3(dna):
    """Takes a DNA sequence and search for a stop codon 'TAG' or 'TAA' or 'TGA' in multiples of 3.
    If there is a stop codon, returns index of first letter
    of the stop codon. If there is no stop codon, return the index of the last letter in the sequence
    >>> find_stop_codon_3('ATTTAGCC')
    3
    >>> find_stop_codon_3('ATTCTATAAGCC')
    6
    >>> find_stop_codon_3('CTACTGTGATGACC')
    6
    >>> find_stop_codon_3('CCTCGGCTACG')
    11
    """
    index = 0
    while index < len(dna):
        if dna[index:index+3] == 'TAG':
            return index
        elif dna[index:index+3] == 'TAA':
            return index
        elif dna[index:index+3] == 'TGA':
            return index
        index = index + 3
    return len(dna)

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATTGACGATGATAGATGGATGTGCCCTAAGCC")
    ['ATGCAT', 'ATGGATGTGCCC']
    >>> find_all_ORFs_oneframe("CCCAGGTTGATGCATTGACGATGATAGATGGATGTGCCCTAAGCC")
    ['ATGCAT', 'ATGGATGTGCCC']
    """
    # TODO: impldna[index_all::]ement this
    index_all = 0
    dna_list = []
    while index_all < len(dna):
        #print(index_all)
        #print(dna[index_all::])
        dna_search = dna[index_all:]
        index = find_start_codon_3(dna_search)
        if index != None:
            new_dna = dna_search[index:]
            stop_index = find_stop_codon_3(new_dna)
            #print (stop_index)
            #print(stop_index)
            dna_list.append(new_dna[:stop_index])

        else:
            return dna_list
        index_all += stop_index +3 + index
        #print(index_all)
    return dna_list


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """

    # TODO: implement this
    ORFlist = []
    new_list1 = find_all_ORFs_oneframe(dna) #find ORF by reading in multiples of 3
    ORFlist += new_list1
    new_list2 = find_all_ORFs_oneframe(dna[1:]) #find ORF by reading in 1 order after multiples of 3
    ORFlist += new_list2
    new_list3 = find_all_ORFs_oneframe(dna[2:]) #find ORF by reading in 2 order after multiples of 3
    ORFlist += new_list3
    return ORFlist


    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    ORFlist = []
    new_list1 = find_all_ORFs(dna) #read it in 1 direction
    ORFlist += new_list1
    new_dna = get_reverse_complement(dna)
    new_list2 = find_all_ORFs(new_dna)
    ORFlist += new_list2
    return ORFlist
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF('ATCCGAATCTAGCAACAAA')
    '0'
    """
    # TODO: implement this
    ORFlist = find_all_ORFs_both_strands(dna)
    if ORFlist != []:
        #print(list)
        longest_ORF = max(ORFlist, key=len)
    else:
        return '0'
    return longest_ORF
    pass



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    ORFlist = []
    for num in range(num_trials):
        new_dna = shuffle_string(dna)
        #print("new_dna", new_dna)
        ORF_longest = longest_ORF(new_dna)
        #print("ORF_longest", ORF_longest)
        ORFlist.append(ORF_longest)
        #print("list", list)
    if ORFlist != []:
        #print("ORFlist", ORFlist)
        longest_ORF_noncoding = max(ORFlist, key = len)
        #print("LongestORFNoncoding", longest_ORF_noncoding)
    else:
        return 0
    return len(longest_ORF_noncoding)

    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    index = 0
    AA_sequence = ''
    AA_list = []
    while index + 3 <= len(dna):

        amino_acid = aa_table[dna[index:index+3]]
        AA_list. append(amino_acid)
        index += 3
    AA = AA_sequence.join(AA_list)
    return AA

    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    AA_list = []
    index = 0
    threshold = longest_ORF_noncoding(dna, 1500)
    print('threshold', threshold)
    ORFlist = find_all_ORFs_both_strands(dna)
    for i in ORFlist:
        if len(i) > threshold:
            AA = coding_strand_to_AA(i)
            AA_list.append(AA)
    return AA_list
    pass

if __name__ == "__main__":
    import doctest
    from amino_acids import aa, codons, aa_table   # you may find these useful
    from load import load_seq
    #doctest.testmod()
    dna = load_seq("./data/X73525.fa")
    result = gene_finder(dna)
    print(result)
    #print(result)
