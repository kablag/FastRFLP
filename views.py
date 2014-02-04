from django.shortcuts import render
from fastrflp.models import PrototypeEnzyme as Pe
from fastrflp.core import Snp
from django.db.models import Q

def index(request):
    return render(request, 'fastrflp/index.html')

def results(request):

    snp = Snp(request.POST['sequence'])
    suppliers = request.POST['suppliers']
    max_num_of_mismatches = int(request.POST.get('max_mismatches'))
    snp.choose_res_to_digest(max_num_of_mismatches=max_num_of_mismatches,
                             suppliers=suppliers)
    result_pes = []
    class ResultPe():
        def __init__(self, name, positions, re_list):
            self.name = name
            self.positions = positions
            self.re_list = re_list


    for pe in snp.penzymes_to_digest:
        cur_pe = Pe.objects.get(name=pe.pe_name)
        query = ''
        if suppliers == '':
            query = "|Q(suppliers__contains = '')"
        for supplier in suppliers:
            query += "|Q(suppliers__contains = '{}')".format(supplier)
        query = query[1:]
        res_for_this_pe = eval("cur_pe.restrictionenzyme_set.filter({})".format(query))
        result_pes.append(ResultPe(pe.pe_name,
                                   pe.positions,
                                   res_for_this_pe))
    return render(request, 'fastrflp/results.html', {
        #  'pes':snp.penzymes_to_digest,
        'result':result_pes,
        })