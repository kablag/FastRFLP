from django.shortcuts import render
from fastrflp.models import PrototypeEnzyme as Pe
from fastrflp.snp import Snp

from django.db.models import Q

def index(request):
    return render(request, 'fastrflp/index.html')

def results(request):
    snp = Snp(request.POST['sequence'])
    suppliers = request.POST['suppliers']
    max_num_of_mismatches = int(request.POST.get('max_mismatches'))
    # snp.choose_res_to_digest(max_num_of_mismatches=max_num_of_mismatches,
    #                          suppliers=suppliers)
    result_pes = []
    class ResultPe():
        def __init__(self, name, positions,outside_snp, re_list):
            #  class Position():
            #      def __init__(self, site):

            self.name = name
            wt_positions = [];
            for wt_pos in positions[0]:
                site_as_list = wt_pos.mismatches_to_list()
                wt_positions.append(site_as_list)
            mut_positions = [];
            for mut_pos in positions[1]:
                mut_positions.append(mut_pos.mismatches_to_list())
            self.wt_positions = wt_positions
            self.mut_positions = mut_positions
            self.wt_outside_snp = outside_snp[0]
            self.mut_outside_snp = outside_snp[1]
            self.re_list = re_list
    pes = Pe.objects.supplied_by(suppliers).can_determine_snp(snp, max_num_of_mismatches)
    for pe in pes:
        result_pes.append(ResultPe(pe[0].name,pe[1],pe[2],pe[0].restrictionenzyme_set.supplied_by(suppliers)))
    return render(request, 'fastrflp/results.html', {
        'snp':snp,
        'result':result_pes,
        })