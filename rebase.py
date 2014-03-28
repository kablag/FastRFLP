from collections import defaultdict

__author__ = 'kablag'

import fastrflp.exceptions as fexceps
import urllib

from fastrflp.models import RestrictionEnzyme as Re
from fastrflp.models import PrototypeEnzyme as Pe

from fastrflp.seq_tools import clean_re_sequence

from django.db import transaction




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



