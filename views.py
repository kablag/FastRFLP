from django.shortcuts import render
from fastrflp.core import Snp


def index(request):
    return render(request, 'fastrflp/index.html')

def results(request):
    snp = Snp(request.GET['sequence'])
    snp.choose_res_to_digest(max_num_of_mismatches=int(request.GET['max_mismatches']),
                             suppliers=request.GET['suppliers'])
    return render(request, 'fastrflp/results.html', {
        'pes':snp.penzymes_to_digest,
        })