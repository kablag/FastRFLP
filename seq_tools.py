__author__ = 'kablag'

import re

ambiguity_dict = { 'a': frozenset('a'),
                   't': frozenset('t'),
                   'g': frozenset('g'),
                   'c': frozenset('c'),

                   'r': frozenset('ga'),
                   'y': frozenset('ct'),
                   'm': frozenset('ac'),
                   'k': frozenset('gt'),
                   's': frozenset('gc'),
                   'w': frozenset('at'),

                   'b': frozenset('cgt'),
                   'd': frozenset('agt'),
                   'h': frozenset('act'),
                   'v': frozenset('acg'),

                   'n': frozenset('atgc'),
                   }

ambiguity_dict_reverse = { frozenset('a')   :'a',
                           frozenset('t')   :'t',
                           frozenset('g')   :'g',
                           frozenset('c')   :'c',

                           frozenset('ga')  :'r',
                           frozenset('ct')  :'y',
                           frozenset('ac')  :'m',
                           frozenset('gt')  :'k',
                           frozenset('gc')  :'s',
                           frozenset('at')  :'w',

                           frozenset('cgt') :'b',
                           frozenset('agt') :'d',
                           frozenset('act') :'h',
                           frozenset('acg') :'v',

                           frozenset('atgc'):'n',
                           }

def expand_sequence(sequence:str):
    esequence = [ambiguity_dict[nuc] for nuc in sequence]
    return esequence

def clean_re_sequence(sequence: str):
    """
    Removes from sequence any chars except "atgcrymkswbdhvn" (case independent)
    :param sequence: Nucleotide sequence
    :return: Cleaned nucleotide sequence

    Example
    =======
    >> clean_sequence('A!aX')
    'Aa'
    """
    return clean_sequence(sequence, 'atgcrymkswbdhvnATGCRYMKSWBDHVN').lower()


def clean_snp_sequence(sequence: str):
    """
    Removes from sequence any chars except "atgcrymkswbdhvn" (case independent)
    :param sequence: Nucleotide sequence
    :return: Cleaned nucleotide sequence

    Example
    =======
    >> clean_sequence('A!aX')
    'Aa'
    """
    return clean_sequence(sequence, 'atgcrymkswbdhvnATGCRYMKSWBDHVN\[\]\/')


def clean_sequence(sequence: str, filter):
    """
    Removes from sequence any chars except filter (case independent)
    :param sequence: Nucleotide sequence
    :return: Cleaned nucleotide sequence

    Example
    =======
    >> clean_sequence('A!aX')
    'Aa'
    """
    return re.sub('[^{}]'.format(filter), '', sequence)


def complement(sequence: str):
    """
    Reads sequence complement

    Complement table (for upper and lower cases)
    =============================================
    A = T
    T = A
    G = C
    C = G
    R (G or A) = Y
    Y (C or T) = R
    M (A or C) = K
    K (G or T) = M
    S (G or C) = S
    W (A or T) = W
    B (not A (C or G or T)) = V
    D (not C (A or G or T)) = H
    H (not G (A or C or T)) = D
    V (not T (A or C or G)) = B
    N (A or C or G or T) = N

    :param sequence: sequence string
    :return: complement sequence string

    Example
    =======
    >> complement('ATGCatgcRYMKSWBDHVNrymkswbdhvn')
    'TACGtacgYRKMSWVHDBNyrkmswvhdbn'
    """
    table = "".maketrans("ATGCatgcRYMKSWBDHVNrymkswbdhvn", "TACGtacgYRKMSWVHDBNyrkmswvhdbn")
    return sequence.translate(table)


def reverse_complement(sequence: str):
    """
    Reads sequence reverse complement

    Complement table (for upper and lower cases)
    =============================================
    (read doc for complement function)

    :param sequence: sequence string
    :return: reverse complement sequence string

    Example
    =======
    >> reverse_complement('ATGCatgcRYMKSWBDHVNrymkswbdhvn')
    'nbdhvwsmkryNBDHVWSMKRYgcatGCAT'
    """
    return complement(sequence)[::-1]


