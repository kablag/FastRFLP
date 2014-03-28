__author__ = 'kik'
import fastrflp.exceptions as fexceps
from fastrflp.seq_tools import ambiguity_dict_reverse

class Site():
    def __init__(self, position, mismatches):
        self.position = position
        self.mismatches = mismatches
        # self.mismatches_list = self.mismatches_to_list()


    def __repr__(self):
        return "{} {}".format(self.position, self.mismatches)

    def add_mismatch(self, position:int, nuc:frozenset):
        self.mismatches.append((position, nuc))

    def mismatches_to_list(self):
        mism = []
        for mis in self.mismatches:
            # nucs = []
            # for n in mis[1]:
            #     nucs.append(n)
            mism.append("{0}:{1}".format(mis[0], ambiguity_dict_reverse[mis[1]]))
        return Site(self.position, mism)


def test_intercept(site:Site, other_site:Site):
    if not site.mismatches and not other_site.mismatches:
        raise fexceps.SitesCollisionError
    misses = []
    site_mis_pos = {p[0] for p in site.mismatches}
    other_site_mis_pos = {p[0] for p in other_site.mismatches}
    if not site_mis_pos.issuperset(other_site_mis_pos):
        return site
    for this_mis in site.mismatches:
        dif = this_mis[1]
        for other_mis in other_site.mismatches:
            if this_mis[0] == other_mis[0]:
                dif = dif.difference(other_mis[1])
        if dif:
            misses.append((this_mis[0], dif))
    if len(misses) != len(site.mismatches):
        raise fexceps.SitesCollisionError
    return Site(site.position, misses)

