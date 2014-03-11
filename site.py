__author__ = 'kik'
import fastrflp.exceptions as fexceps


class Site():
    def __init__(self, position, mismatches):
        self.position = position
        self.mismatches = mismatches

    def __repr__(self):
        return "{} {}".format(self.position, self.mismatches)

    def add_mismatch(self, position:int, nuc:frozenset):
        self.mismatches.append((position, nuc))


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

