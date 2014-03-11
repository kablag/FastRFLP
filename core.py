from collections import defaultdict

__author__ = 'kablag'

import fastrflp.exceptions as fexceps
import urllib
import re
import uuid
import time
from fastrflp.models import RestrictionEnzyme as Re
from fastrflp.models import PrototypeEnzyme as Pe
from fastrflp.snp import Snp
from fastrflp.seq_tools import expand_sequence, clean_re_sequence
from django.db.models import Q
from django.db import transaction

def can_recognize(a:frozenset, b:frozenset):
    return not a.isdisjoint(b)

def gen_mask(pe_len, pe_seq, template_seq, snp_pos, allele_pos_end, allele_len,
             max_num_of_mismatches):
    sites = []
    for pe_pos in range(snp_pos - pe_len, allele_pos_end):
        try:
            cur_num_mismatches = 0
            mismatches = []
            for tnuc_pos in range(pe_pos, pe_pos + pe_len):
                if not can_recognize(pe_seq[tnuc_pos - pe_pos],
                                     template_seq[tnuc_pos]):
                    if len(mismatches) < max_num_of_mismatches:
                        mismatch_position = tnuc_pos if tnuc_pos < snp_pos \
                            else tnuc_pos - allele_len
                        mismatches.append((mismatch_position,
                                           pe_seq[tnuc_pos - pe_pos]))
                    else: raise fexceps.DigestError
            sites.append((pe_pos,mismatches))
        except fexceps.DigestError:
            pass
    return sites


def pe_can_determine_snp(clean_recognition_sequence, snp:Snp, max_num_of_mismatches):
    pe_len = len(clean_recognition_sequence)
    pe_seq = expand_sequence(clean_recognition_sequence.lower())
    wt_sites = gen_mask(pe_len,
                        pe_seq,
                        snp.wt_sequence,
                        snp.snp_pos,
                        snp.wt_pos_end,
                        snp.wt_allele_len,
                        max_num_of_mismatches,
                        )
    mut_sites = gen_mask(pe_len,
                        pe_seq,
                        snp.mut_sequence,
                        snp.snp_pos,
                        snp.mut_pos_end,
                        snp.mut_allele_len,
                        max_num_of_mismatches,
                        )
    return (wt_sites, mut_sites)



def can_recognize(a:frozenset, b:frozenset):
    return not a.isdisjoint(b)

def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print('{} ({}, {}) {} sec'.format(method.__name__, args, kw, te - ts))
        return result

    return timed


def nuc_can_recognize(recognition_nucleotide: str, template_nucleotide: str):
    """
            Tests if one nucleotide can recognize another (case independent)

            Only one nucleotide can be tested!!! This function is not for sequences!!!

            :param recognition_nucleotide: Recognition nucleotide
            :param template_nucleotide: Template nucleotide
            :return: True if can recognize

            Example
            =======
            >> nuc_can_recognize('r','a')
            True
            >> nuc_can_recognize('r','c')
            False

            """
    unrecognition_dict = {'a': 'tgc yksb',
                          't': 'agc rmsv',
                          'g': 'atc ymwh',
                          'c': 'atg rkwd',
                          'r': 'ct y',
                          'y': 'ga r',
                          'm': 'tg k',
                          'k': 'ac m',
                          's': 'at w',
                          'w': 'gc s',
                          'b': 'a',
                          'd': 'c',
                          'h': 'g',
                          'v': 't',
                          'n': '', }

    if not len(recognition_nucleotide) == 1:
        raise fexceps.TemplateRecognitionError('Not one nucleotide provided: %r' % recognition_nucleotide)

    recognition_nucleotide = recognition_nucleotide.lower()
    template_nucleotide = template_nucleotide.lower()
    try:
        if recognition_nucleotide not in unrecognition_dict[template_nucleotide]:
            return True
        else:
            return False
    except KeyError:
        raise fexceps.TemplateRecognitionError('No such nucleotide in dict: %r' % template_nucleotide)


def gen_mask(pe_len, pe_seq, template_seq, snp_pos, allele_pos_end, allele_len,
             max_num_of_mismatches):
    sites = []
    for pe_pos in range(snp_pos - pe_len, allele_pos_end):
        try:
            cur_num_mismatches = 0
            mismatches = []
            for tnuc_pos in range(pe_pos, pe_pos + pe_len):
                if not can_recognize(pe_seq[tnuc_pos - pe_pos],
                                     template_seq[tnuc_pos]):
                    if len(mismatches) < max_num_of_mismatches:
                        mismatch_position = tnuc_pos if tnuc_pos < snp_pos \
                            else tnuc_pos - allele_len
                        mismatches.append((mismatch_position,
                                           pe_seq[tnuc_pos - pe_pos]))
                    else: raise fexceps.DigestError
            sites.append((pe_pos,mismatches))
        except fexceps.DigestError:
            pass
    return sites


def pe_can_determine_snp(pe:Pe, snp:Snp, max_num_of_mismatches):
    def get_recognition_mask(re_sequence,
                                 target_sequence,
                                 snp_position,
                                 detect_allele,
                                 opposite_allele):
        """
         :param re_sequence:
        :param target_sequence:
        :param snp_position:
        :param detect_allele:
        :param opposite_allele:
        :return: :raise DigestError:
        """
        if not nuc_can_recognize(re_sequence[snp_position], detect_allele):
            raise fexceps.DigestError
        if nuc_can_recognize(re_sequence[snp_position], opposite_allele):
            raise fexceps.DigestError
        mask = []
        num_of_mismatches = 0
        mismatch_before_snp = False
        for i in range(len(re_sequence)):
            if not nuc_can_recognize(re_sequence[i], target_sequence[i]):
                mask.append((i, re_sequence[i]))
                num_of_mismatches += 1
                if num_of_mismatches > max_num_of_mismatches \
                        or (mismatch_before_snp and i > snp_position) \
                        or i == snp_position:
                    raise fexceps.DigestError
                if i < snp_position:
                    mismatch_before_snp = True
        return mask


    pe_len = len(pe.clean_recognition_sequence)
    pe_seq = expand_sequence(pe.clean_recognition_sequence.lower())
    wt_sites = gen_mask(pe_len,
                        pe_seq,
                        snp.wt_sequence,
                        snp.snp_pos,
                        snp.wt_pos_end,
                        snp.wt_allele_len,
                        max_num_of_mismatches,
                        )
    mut_sites = gen_mask(pe_len,
                        pe_seq,
                        snp.mut_sequence,
                        snp.snp_pos,
                        snp.mut_pos_end,
                        snp.mut_allele_len,
                        max_num_of_mismatches,
                        )
    return (wt_sites, mut_sites)

# class Snp_old():
#     #@timeit
#     @classmethod
#     def choose_res_to_digest(self, max_num_of_mismatches=0, suppliers=''):
#
#
#         def get_recognition_mask(re_sequence,
#                                  target_sequence,
#                                  snp_pos,
#                                  detect_allele,
#                                  opposite_allele):
#             """
#
#             :param re_sequence:
#             :param target_sequence:
#             :param snp_pos:
#             :param detect_allele:
#             :param opposite_allele:
#             :return: :raise DigestError:
#             """
#             if not nuc_can_recognize(re_sequence[snp_pos], detect_allele):
#                 raise fexceps.DigestError
#             if nuc_can_recognize(re_sequence[snp_pos], opposite_allele):
#                 raise fexceps.DigestError
#             mask = []
#             num_of_mismatches = 0
#             mismatch_before_snp = False
#             for i in range(len(re_sequence)):
#                 if not nuc_can_recognize(re_sequence[i], target_sequence[i]):
#                     mask.append((i, re_sequence[i]))
#                     num_of_mismatches += 1
#                     if num_of_mismatches > max_num_of_mismatches \
#                             or (mismatch_before_snp and i > snp_pos) \
#                             or i == snp_pos:
#                         raise fexceps.DigestError
#                     if i < snp_pos:
#                         mismatch_before_snp = True
#             return mask
#
#         class PeToDigest():
#             def __init__(self, pe_name):
#                 self.pe_name = pe_name
#                 self.positions = []
#
#             def add_position(self, position, mask, wt):
#                 self.positions.append((position, mask, wt))
#
#             def __repr__(self):
#                 return '%r' % self.pe_name
#
#         query = ''
#         if suppliers == '':
#             query = "|Q(restrictionenzyme__suppliers__contains = '')"
#         for supplier in suppliers:
#             query += "|Q(restrictionenzyme__suppliers__contains = '{}')".format(supplier)
#         query = query[1:]
#         for penzyme in eval("Pe.objects.filter({}).distinct()".format(query)):
#             cur_pe = PeToDigest(penzyme.name)
#             l = len(penzyme.clean_recognition_sequence)
#
#             for i in range(self.snp_pos - l,
#                            self.snp_pos):
#                 try:
#                     mask = get_recognition_mask(penzyme.clean_recognition_sequence,
#                                                 self.wt_sequence[i:i + l],
#                                                 self.snp_pos - i - 1,
#                                                 self.wt_allele,
#                                                 self.mut_allele)
#                     cur_pe.add_position(i, mask, True)
#                 except fexceps.DigestError:
#                     pass
#                 try:
#                     mask = get_recognition_mask(penzyme.clean_recognition_sequence,
#                                                 self.mut_sequence[i:i + l],
#                                                 self.snp_pos - i - 1,
#                                                 self.mut_allele,
#                                                 self.wt_allele)
#                     cur_pe.add_position(i, mask, False)
#                 except fexceps.DigestError:
#                     pass
#
#             if cur_pe.positions:
#                 self.penzymes_to_digest.append(cur_pe)
#
#     def __init__(self, sequence: str, name=None):
#         def from_IUPAC(sequence: str):
#             snp_positions = []
#             for char in "rymksw":
#                 if sequence.count(char) == 1:
#                     position = sequence.find(char)
#                     snp_positions.append((char, position, IUPACdict[char]))
#             if len(snp_positions) == 1:
#                 return snp_positions[0]
#             return None
#
#         def from_wt_slash_mut(sequence: str):
#             snps = re.findall('\[[atgc]/[atgc]]', sequence)
#             if len(snps) == 1:
#                 return snps[0], sequence.find(snps[0]), \
#                        tuple(re.findall('[atgc]', snps[0]))
#             return None
#
#         self.name = name if name is not None else str(uuid.uuid1())
#         sequence = clean_snp_sequence(sequence.lower())
#         info_from_seq = from_wt_slash_mut(sequence)
#         if info_from_seq is None:
#             info_from_seq = from_IUPAC(sequence)
#         if info_from_seq is not None:
#             self.original_snp_sign, self.snp_pos, \
#             (self.wt_allele, self.mut_allele) = info_from_seq
#             self.snp_pos += 1  # from 0 based to real positions
#             self.wt_sequence = sequence.replace(self.original_snp_sign, self.wt_allele, 1)
#             self.mut_sequence = sequence.replace(self.original_snp_sign, self.mut_allele, 1)
#             self.digest_penzymes = []
#         else:
#             raise fexceps.GetSNPFromSequenceError('No valid SNP info in seq: %r' % sequence)
#
#     def __repr__(self):
#         return self.name
#
#     def __str__(self):
#         return (
#             """
# {name}
# SNP postition: {snp_pos}
# WT: {wt_allele}
# {wt_seq}
# MUT: {mut_allele}
# {mut_seq}
# """.format(name=self.name,
#            snp_pos=self.snp_pos,
#            wt_allele=self.wt_allele,
#            wt_seq=self.wt_sequence,
#            mut_allele=self.mut_allele,
#            mut_seq=self.mut_sequence, ))
#

@transaction.atomic
def update_rebase_from_url(url='http://rebase.neb.com/rebase/link_itype2'):
    """
    Updates (or creates) REs database from url (NEB Type II format with tabs)

    Table Format
    ============
    enzyme name [tab] prototype [tab] recognition sequence with cleavage site
    [tab] methylation site and type [tab] commercial source [tab] references

    :param url: REs table location (default is 'http://rebase.neb.com/rebase/link_itype2')
    :returns: list of new REs and list of updated REs

    Example
    =======
    >> list_of_new_res, list_of_updated_res = update_rebase_from_url()
    Downloaded REs Entries: 3725
    New REs: 3725
    Updated REs: 0
    """
    # from fastrflp.models import PrototypeEnzyme as Pe
    # from fastrflp.models import RestrictionEnzyme as Re

    

    try:
        fh = urllib.request.urlopen(url)
    except (ValueError, urllib.error.URLError):
        raise fexceps.UpdateRebaseError('Error in URL: %r' % url)
    html = fh.read().decode('utf8')
    new_pes = []
    upd_pes = []
    new_res = []
    upd_res = []
    try:
        enzyme_rows = html.splitlines()[10:-1]
    except IndexError:
        raise fexceps.UpdateRebaseError('Document %r is to short' % url)
    for enzyme_row in enzyme_rows:
        # enzyme name [tab] prototype [tab] recognition sequence with cleavage site
        # [tab] methylation site and type [tab] commercial source [tab] references
        try:
            re_name, re_prototype, re_recognition_sequence, re_methylation_site, re_suppliers, re_refs \
                = tuple(enzyme_row.split('\t'))
        except ValueError as err:
            if 'unpack' in err:
                raise fexceps.UpdateRebaseError('RE table format error (num of columns)')
        if re_prototype == '':
            re_prototype = re_name

        if not Pe.objects.filter(name=re_prototype).exists():
            new_pe = Pe(name=re_prototype,
                        clean_recognition_sequence=clean_re_sequence(re_recognition_sequence),
                        )
            new_pe.save()
            new_pe.restrictionenzyme_set.create(name=re_name,
                                                recognition_sequence=re_recognition_sequence,
                                                suppliers=re_suppliers)
            new_res.append(re_name)
            #new_pe.save()
            new_pes.append(re_prototype)
        else:
            # update PE
            existing_pe = Pe.objects.get(name=re_prototype)
            if existing_pe.clean_recognition_sequence != clean_re_sequence(re_recognition_sequence):
                # existing_pe.prototype = re_prototype
                # existing_pe.recognition_sequence = re_recognition_sequence
                existing_pe.clean_recognition_sequence = clean_re_sequence(re_recognition_sequence)
                # existing_pe.suppliers = re_suppliers
                existing_pe.save()
                upd_pes.append(re_prototype)
            if not existing_pe.restrictionenzyme_set.filter(name=re_name).exists():
                existing_pe.restrictionenzyme_set.create(name=re_name,
                                                         recognition_sequence=re_recognition_sequence,
                                                         suppliers=re_suppliers)
                #  existing_pe.save()
                new_res.append(re_name)
            else:
                existing_re = Re.objects.get(name=re_name)
                if existing_re.recognition_sequence != re_recognition_sequence \
                or existing_re.suppliers != re_suppliers:
                    existing_re.recognition_sequence = re_recognition_sequence
                    existing_re.suppliers = re_suppliers
                    existing_re.save()
                    upd_res.append(re_name)
    print('''
          Downloaded REs Entries: {}
          New PEs: {}\tUpdated PEs: {}
          New REs: {}\tUpdated REs: {}
          '''.format(len(enzyme_rows),
                     len(new_pes), len(upd_pes),
                     len(new_res), len(upd_res)))
    return new_pes, upd_pes, new_res, upd_res



