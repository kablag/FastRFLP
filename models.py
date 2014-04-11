from django.db import models
#from django.db.models.query import QuerySet
from django.db.models import Q
from fastrflp.snp import Snp
import fastrflp.exceptions as fexceps
from fastrflp.seq_tools import expand_sequence
from fastrflp.seq_tools import reverse_complement
from fastrflp.site import Site
from fastrflp.site import test_intercept
from fastrflp.timeit import timeit

#  class PrototypeEnzymeManager(models.Manager):

class QuerySetManager(models.Manager):
    def get_query_set(self):
        model = models.get_model(self.model._meta.app_label, self.model._meta.object_name)
        return PeQuerySet(model)

    def __getattr__(self, attr, *args):
        try:
            return getattr(self.__class__, attr, *args)
        except AttributeError:
            return getattr(self.get_query_set(), attr, *args)

class PeQuerySet(models.query.QuerySet):
    def can_determine_snp(self, snp:Snp, max_num_of_mismatches, do_not_allow_mismatches_near_snp):
        def can_recognize(a:frozenset, b:frozenset):
            return not a.isdisjoint(b)

        def gen_mask(pe_len, pe_seq, template_seq, snp_pos,
                     allele_pos_end, allele_len):
            sites = []
            exclude_poss = range(snp_pos - 2, allele_pos_end + 1) if do_not_allow_mismatches_near_snp \
                else range(snp_pos - 1, allele_pos_end)
            for pe_pos in range(snp_pos - pe_len, allele_pos_end):
                try:
                    mismatches = set()
                    for tnuc_pos in range(pe_pos, pe_pos + pe_len):
                        if not can_recognize(pe_seq[tnuc_pos - pe_pos],
                                             template_seq[tnuc_pos]):
                            if len(mismatches) < max_num_of_mismatches and \
                                            tnuc_pos not in exclude_poss:
                                mismatch_position = tnuc_pos if tnuc_pos < snp_pos \
                                    else tnuc_pos - allele_len
                                mismatches.add((mismatch_position,
                                                   pe_seq[tnuc_pos - pe_pos]))
                            else:
                                raise fexceps.DigestError
                    sites.append(Site(pe_pos,mismatches))
                except fexceps.DigestError:
                    pass
            return sites

        def pe_can_determine_snp(pe:PrototypeEnzyme,
                                 snp:Snp,
                                 max_num_of_mismatches):
            pe_len = len(pe.clean_recognition_sequence)
            pe_seq = expand_sequence(pe.clean_recognition_sequence.lower())
            wt_sites = gen_mask(pe_len,
                                pe_seq,
                                snp.ex_wt_sequence,
                                snp.snp_pos,
                                snp.wt_pos_end,
                                snp.wt_allele_len,
                                )
            mut_sites = gen_mask(pe_len,
                                 pe_seq,
                                 snp.ex_mut_sequence,
                                 snp.snp_pos,
                                 snp.mut_pos_end,
                                 snp.mut_allele_len,
                                 )
            if not pe.clean_recognition_sequence == \
                reverse_complement(pe.clean_recognition_sequence):
                pe_seq_rc = expand_sequence(reverse_complement(pe.clean_recognition_sequence.lower()))
                wt_sites = wt_sites + gen_mask(pe_len,
                                pe_seq_rc,
                                snp.ex_wt_sequence,
                                snp.snp_pos,
                                snp.wt_pos_end,
                                snp.wt_allele_len,
                                )
                mut_sites = mut_sites + gen_mask(pe_len,
                                 pe_seq_rc,
                                 snp.ex_mut_sequence,
                                 snp.snp_pos,
                                 snp.mut_pos_end,
                                 snp.mut_allele_len,
                                 )
            wt_filtered = []
            for wt_site in wt_sites:
                try:
                    for mut_site in mut_sites:
                        wt_site = test_intercept(wt_site, mut_site)
                    wt_filtered.append(wt_site)
                except fexceps.SitesCollisionError:
                    pass

            mut_filtered = []
            for mut_site in mut_sites:
                try:
                    for wt_site in wt_sites:
                        mut_site = test_intercept(mut_site, wt_site)
                    mut_filtered.append(mut_site)
                except fexceps.SitesCollisionError:
                    pass

            return (wt_filtered, mut_filtered)

        d_list = []
        for pe in self.all():
            can_determine = pe_can_determine_snp(pe,
                                                 snp,
                                                 max_num_of_mismatches)
            if can_determine[0] or can_determine[1]:
                def sites_outside_snp(target_seq, snp_pos, allele_pos_end):
                    sites_outside = []
                    pe_len = len(pe.clean_recognition_sequence)
                    pe_seq = expand_sequence(pe.clean_recognition_sequence.lower())
                    # seq_before_snp = snp.ex_wt_sequence[:snp.snp_pos - pe_len]
                    # seq_after_snp = snp.ex_wt_sequence[snp.snp_pos + snp.wt_allele_len:]
                    for i in range(0, len(target_seq) - pe_len):
                        try:
                            if i in range(snp_pos - pe_len, allele_pos_end):
                                raise fexceps.DigestError
                            for ii in range(0, pe_len):
                                if not can_recognize(target_seq[i+ii], pe_seq[ii]):
                                    raise fexceps.DigestError
                            sites_outside.append(i)
                        except fexceps.DigestError:
                            pass
                    if not pe.clean_recognition_sequence == \
                        reverse_complement(pe.clean_recognition_sequence):
                        pe_seq_rc = expand_sequence(reverse_complement(pe.clean_recognition_sequence.lower()))
                        for i in range(0, len(target_seq) - pe_len):
                            try:
                                if i in range(snp_pos - pe_len, allele_pos_end):
                                    raise fexceps.DigestError
                                for ii in range(0, pe_len):
                                    if not can_recognize(target_seq[i+ii], pe_seq_rc[ii]):
                                        raise fexceps.DigestError
                                sites_outside.append(i)
                            except fexceps.DigestError:
                                pass
                    return sites_outside
                wt_sites_outside_snp = sites_outside_snp(snp.ex_wt_sequence, snp.snp_pos, snp.wt_pos_end)
                mut_sites_outside_snp = sites_outside_snp(snp.ex_mut_sequence, snp.snp_pos, snp.mut_pos_end)

                d_list.append((pe,
                               can_determine,
                               (wt_sites_outside_snp,
                                mut_sites_outside_snp)))
        return d_list

    def supplied_by(self, suppliers:str):
        query_elms = []
        if suppliers == '':
            query_elms.append("Q(restrictionenzyme__suppliers__contains = '')")
            query_elms.append("|")
        for supplier in suppliers:
            query_elms.append("Q(restrictionenzyme__suppliers__contains = '{}')"
                              .format(supplier))
            query_elms.append("|")
        query = "".join(query_elms[:-1])
        return eval("self.filter({})".format(query)).distinct()


class PrototypeEnzyme(models.Model):
    name = models.CharField(max_length=15)
    clean_recognition_sequence = models.CharField(max_length=30)
     #  object = PrototypeEnzymeManager()
    objects = QuerySetManager()
    def __str__(self):
        res = []
        for re in self.restrictionenzyme_set.all():
            res.append('\t{} {} {}\n'.format(re.name, re.recognition_sequence, re.suppliers))
        return '{} {} \n{}'\
            .format(self.name,
                    self.clean_recognition_sequence,
                    ''.join(res)) if res \
            else '{} {}'\
            .format(self.name,
                    self.clean_recognition_sequence)

class RestrictionEnzymeManager(models.Manager):
    def supplied_by(self, suppliers:str):
        query_elms = []
        if suppliers == '':
            query_elms.append("Q(suppliers__contains = '')")
            query_elms.append("|")
        for supplier in suppliers:
            query_elms.append("Q(suppliers__contains = '{}')".format(supplier))
            query_elms.append("|")
        query = "".join(query_elms[:-1])
        return eval("self.filter({})".format(query))

class RestrictionEnzyme(models.Model):
    prototype = models.ForeignKey(PrototypeEnzyme)
    name = models.CharField(max_length=15)
    recognition_sequence = models.CharField(max_length=30)
    suppliers = models.CharField(max_length=20)
    objects = RestrictionEnzymeManager()

    def __str__(self):
        return '{} {} {}'.format(self.name, self.recognition_sequence, self.suppliers)



