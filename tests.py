from django.test import TestCase

from fastrflp.models import PrototypeEnzyme as Pe

from fastrflp.snp import Snp
import fastrflp.exceptions as fexceps
from fastrflp.seq_tools import expand_sequence

# Create your tests here.


class SnpClassTests(TestCase):
    def test_init_from_IUPAC(self):
        """
        Snp("aaaaaaaaaaaaaaaaaaaaaaYaaaaaaaaaaaaaaaaaaaaaaa" ok
        """
        snp = Snp("aaaaaaaaaaaaaaaaaaaaaaYaaaaaaaaaaaaaaaaaaaaaaa")
        self.assertEqual(snp.snp_pos, 23)
        self.assertEqual(snp.wt_allele, 'c')
        self.assertEqual(snp.ex_wt_sequence,
                         expand_sequence('aaaaaaaaaaaaaaaaaaaaaacaaaaaaaaaaaaaaaaaaaaaaa'))
        self.assertEqual(snp.mut_allele, 't')
        self.assertEqual(snp.ex_mut_sequence,
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