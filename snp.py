__author__ = 'kablag'
import re
import uuid
import fastrflp.exceptions as fexceps
from fastrflp.seq_tools import expand_sequence
from fastrflp.seq_tools import ambiguity_dict

class Snp():
    def __init__(self, sequence: str, name=None):
        def from_IUPAC(sequence: str):
            snp_positions = []
            for char in "rymksw":
                if sequence.count(char) == 1:
                    position = sequence.find(char)
                    snp_positions.append((char, position, tuple(ambiguity_dict[char])))#  IUPACdict[char]))
            if len(snp_positions) == 1:
                return snp_positions[0]
            return None

        def from_wt_slash_mut(sequence: str):
            snps = re.findall('\[[atgc]*/[atgc]*]', sequence)
            if len(snps) == 1:
                return snps[0], sequence.find(snps[0]), \
                       tuple(snps[0]
                             .replace('[','')
                             .replace(']','')
                             .split('/'))
            return None

        self.name = name if name is not None else str(uuid.uuid1())
        sequence = re.sub('[^{}]'.format('atgcrymkswbdhvnATGCRYMKSWBDHVN\[\]\/'),
                          '',
                          sequence.lower())
        info_from_seq = from_wt_slash_mut(sequence)
        if info_from_seq is None:
            info_from_seq = from_IUPAC(sequence)
        if info_from_seq is not None:
            self.original_snp_sign, self.snp_pos, \
            (self.wt_allele, self.mut_allele) = info_from_seq
            self.snp_pos += 1  # from 0 based to real positions
            self.wt_allele_len = len(self.wt_allele)
            self.wt_pos_end = self.snp_pos + self.wt_allele_len  - 1
            self.mut_allele_len = len(self.mut_allele)
            self.mut_pos_end = self.snp_pos + self.mut_allele_len - 1
            self.wt_sequence = expand_sequence(sequence.replace(self.original_snp_sign,
                                                                self.wt_allele, 1))
            self.mut_sequence = expand_sequence(sequence.replace(self.original_snp_sign,
                                                                 self.mut_allele, 1))
            self.digest_penzymes = []
        else:
            raise fexceps.GetSNPFromSequenceError('No valid SNP info in seq: %r' % sequence)
    def __repr__(self):
        return self.name

    def __str__(self):
        return (
            """
{name}
SNP postition: {snp_pos}
WT: {wt_allele} wt allele pos end:{wt_pos_end}
{wt_seq}
MUT: {mut_allele} mut allele pos end:{mut_pos_end}
{mut_seq}
""".format(name=self.name,
           snp_pos=self.snp_pos,
           wt_allele=self.wt_allele,
           wt_pos_end=self.wt_pos_end,
           wt_seq=self.wt_sequence,
           mut_allele=self.mut_allele,
           mut_pos_end=self.mut_pos_end,
           mut_seq=self.mut_sequence, ))