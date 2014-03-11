from django.test import TestCase
import fastrflp

from fastrflp.models import PrototypeEnzyme as Pe
from fastrflp.core import nuc_can_recognize
#from fastrflp.core import Snp
from fastrflp.core import Snp
import fastrflp.exceptions as fexceps
from fastrflp.core import pe_can_determine_snp as pd
from fastrflp.seq_tools import expand_sequence
import string
import re
# Create your tests here.


class NucCanRecognizeTests(TestCase):
    def test_recognition(self):
        """
        nuc_can_recognize() should return False for all cases with
        unrecognition_dictionary
        """
        unrecognition_dict = {'a': 'tgcyksb',
                                  't': 'agcrmsv',
                                  'g': 'atcymwh',
                                  'c': 'atgrkwd',
                                  'r': 'cty',
                                  'y': 'gar',
                                  'm': 'tgk',
                                  'k': 'acm',
                                  's': 'atw',
                                  'w': 'gcs',
                                  'b': 'a',
                                  'd': 'c',
                                  'h': 'g',
                                  'v': 't',
                                  'n': '', }
        ok = True
        for recognition_nucleotide in unrecognition_dict:
            for template_nucleotide in unrecognition_dict[recognition_nucleotide]:
                if nuc_can_recognize(recognition_nucleotide, template_nucleotide):
                    self.assertFalse(False, '"{}" can recognize "{}"'.format(recognition_nucleotide,
                                                                         template_nucleotide))

    def test_template_nucleotide_not_in_dict(self):
        """
        nuc_can_recognize(ANY, X), where X not in ATGCatgcRYMKSWBDHVNrymkswbdhvn
        should return TemplateRecognitionError('No such nucleotide in dict')
        """
        for c in re.sub('[ATGCatgcRYMKSWBDHVNrymkswbdhvn]', '', string.printable):
            try:
                self.assertTrue(nuc_can_recognize(' ', c), '{} in recognition dictionary'.format(c))
            except fexceps.TemplateRecognitionError as err:
                self.assertFalse('No such nucleotide in dict' not in err.code)

    def test_no_recognition_nucleotides_provided(self):
        """
        nuc_can_recognize(X, ANY), when len(X) != 1
        should return TemplateRecognitionError('Not one nucleotide provided')
        """
        try:
            nuc_can_recognize('','a')
        except fexceps.TemplateRecognitionError as err:
            self.assertFalse('Not one nucleotide provided' not in err.code)

    def test_more_than_one_recognition_nucleotide_provided(self):
        """
        nuc_can_recognize(X, ANY), when len(X) != 1
        should return TemplateRecognitionError('Not one nucleotide provided')
        """
        try:
            nuc_can_recognize('  ','a')
        except fexceps.TemplateRecognitionError as err:
            self.assertFalse('Not one nucleotide provided' not in err.code)


class SnpClassTests(TestCase):
    def test_init_from_IUPAC(self):
        """
        Snp("aaaaaaaaaaaaaaaaaaaaaaYaaaaaaaaaaaaaaaaaaaaaaa" ok
        """
        snp = Snp("aaaaaaaaaaaaaaaaaaaaaaYaaaaaaaaaaaaaaaaaaaaaaa")
        self.assertEqual(snp.snp_pos, 23)
        self.assertEqual(snp.wt_allele, 'c')
        self.assertEqual(snp.wt_sequence,
                         expand_sequence('aaaaaaaaaaaaaaaaaaaaaacaaaaaaaaaaaaaaaaaaaaaaa'))
        self.assertEqual(snp.mut_allele, 't')
        self.assertEqual(snp.mut_sequence,
                         expand_sequence('aaaaaaaaaaaaaaaaaaaaaataaaaaaaaaaaaaaaaaaaaaaa'))

    def test_init_from_wt_slash_mut(self):
        """
        Snp("aaaaaaaaaaaaaaaaaaaaaa[a/g]aaaaaaaaaaaaaaaaaaaaaaa" ok
        """
        Snp("aaaaaaaaaaaaaaaaaaaaaa[a/g]aaaaaaaaaaaaaaaaaaaaaaa")

    def test_init_from_IUPAC_with_two_snps_in_seq(self):
        """
        Snp("aaaaaaaaaaaaaaaaaaaaaaYaaaaaaYaaaaaaaaaaaaaaaa" Error
        """
        try:
            Snp("aaaaaaaaaaaaaaaaaaaaaaYaaaaaaYaaaaaaaaaaaaaaaa")
        except fexceps.GetSNPFromSequenceError as err:
            self.assertTrue('No valid SNP info in seq' in err.code)

class DetermitionTests(TestCase):
    def test_determine_without_mismatches(self):
        p = Pe()
        p.name = 'TaqI'
        p.clean_recognition_sequence = 'tcga'
        p.save()
        s = Snp('aaaaaat[tcgatc/c]gaaaaa')
        print(Pe.objects.can_determine_snp(s, 0))

    def test_determine_when_mismatches(self):
        p = Pe()
        p.name = 'TaqI'
        p.clean_recognition_sequence = 'tcga'
        p.save()
        s = Snp('aaaaaat[tcgatc/c]gaaaaa')
        print(Pe.objects.can_determine_snp(s, 0))

    def test_determine_without_mismatches_IUPAC(self):
        p = Pe()
        p.name = 'TaqI'
        p.clean_recognition_sequence = 'tcga'
        p.save()
        s = Snp('aaaaaatycgaaaaa')
        print("{} \n{}".format(self._testMethodName,
                               Pe.objects.can_determine_snp(s, 0)))

    def test_determine_with_one_mismatch(self):
        p = Pe()
        p.name = 'TaqI'
        p.clean_recognition_sequence = 'tcga'
        p.save()
        s = Snp('aaaaaaa[tc/c]gaaaaa')
        print(Pe.objects.can_determine_snp(s, 1))

    def test_determine_with_two_mismatches(self):
        p = Pe()
        p.name = 'TaqI'
        p.clean_recognition_sequence = 'tcga'
        p.save()
        s = Snp('aaaaaaa[tc/c]aaaaaa')
        print(Pe.objects.all().can_determine_snp(s, 2))

    def test_not_determine_with_mismatches_in_snp(self):
        p = Pe()
        p.name = 'TaqI'
        p.clean_recognition_sequence = 'tcga'
        p.save()
        s = Snp('aaaaaat[g/t]gaaaaa')
        print(Pe.objects.all().can_determine_snp(s, 1))

    def test_not_determine_when_wt_mis_generate_mut_site(self):
        p = Pe()
        p.name = 'TaqI'
        p.clean_recognition_sequence = 'tcga'
        p.save()
        s = Snp('aaaaaat[g/t]gaaaaa')
        print(Pe.objects.all().can_determine_snp(s, 1))

    def test_3(self):
        p = Pe()
        p.name = 'TaqI'
        p.clean_recognition_sequence = 'tccy'
        p.save()
        s = Snp('aaaaaat[c/t]cataaa')
        print("{} \n{}".format(self._testMethodName,
                               Pe.objects.can_determine_snp(s, 1)))