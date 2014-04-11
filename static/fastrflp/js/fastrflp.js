/**
 * Created by kablag on 11.03.14.
 */
function cleanseqlen(seq)
{
    return seq.replace(/[^atgcrymkswbdhvnATGCRYMKSWBDHVN]/g, "").length;
}

function colorizeseq(snp, e)
{
    e = e || window.event;
    var data = [];
    var target = e.srcElement || e.target;
    while (target && target.nodeName !== "TR") {
        target = target.parentNode;
    }
    if (target) {
        var cells = target.getElementsByTagName("td");
        for (var i = 0; i < cells.length; i++) {
            data.push(cells[i]);
        }
    }
    var enzname = data[1].innerHTML;
    var enzseq = data[2].innerHTML;
    var enzlen = cleanseqlen(enzseq);
    var alleletype = data[4].innerHTML;
    var enzpos = Number(data[5].innerHTML);
    var mismatches = data[6].innerHTML.split('@')[1];
    var outsidesnp = data[7].title.match(/\d+/g);;

    var mask = [];
    var snppos = snp.snp_pos;
    if(alleletype=='Wt')
    {
        var seq = snp.wt_seq;
        var snpposend = snp.wt_pos_end;
        var allelelen = snp.wt_allele_len;
    }
    else
    {
        var seq = snp.mut_seq;
        var snpposend = snp.mut_pos_end;
        var allelelen = snp.mut_allele_len;
    }
    for(var i = 0; i < seq.length; i++)
    {
        mask.push('-');
    };
    // add snp
    for (var i = snppos - 1; i < snpposend; i++) {
        mask[i]= '<span class="snppos">'+ seq[i] + '</span>';
    }
    if(e){
        // add enzime
        for (var i = enzpos; i < enzpos + enzlen ; i++) {
            if(mask[i]=="-") {
                mask[i]='<span class="enzpos">'+ seq[i] + '</span>';
            }
        }

        // add mism
        if(mismatches) {
            var mis_posns = mismatches.match(/\d+/g);
            var mis_nucs = mismatches.match(/[a-z]/g);
            function addmismatch(pos, index)
            {
                if(pos >= snppos - 1) {
                    mask[parseInt(pos)+allelelen]='<span class="misenz"'
                        + 'title="'
                        + seq[parseInt(pos)+allelelen]
                        + '">'
                        + mis_nucs[index]
                        + '</span>';
                }
                else {
                    mask[parseInt(pos)]='<span class="misenz"'
                        + 'title="'
                        + seq[parseInt(pos)]
                        + '">'
                        + mis_nucs[index]
                        + '</span>';
                }
            };
            mis_posns.forEach(addmismatch);
        }
        // add outside snp digest
        if(outsidesnp) {
            function addoutsidesnp(pos) {
                pos = Number(pos);
                for (var i = pos; i < pos + enzlen ; i++) {
                    if(mask[i]=="-") {
                        mask[i]='<span class="outsnp">'+ seq[i] + '</span>';
                    };
                };
            };
            outsidesnp.forEach(addoutsidesnp);
        }
    }
    var resultseq = '';
    for (var i = 0; i < mask.length ; i++) {
        if(mask[i]=='-'){
            resultseq += seq[i];
        }
        else
        {
            resultseq += mask[i];
        }
    }
    document.getElementById("snpsequence").innerHTML = resultseq;
    if(mismatches){
        document.getElementById("selectedre").innerHTML = "Enzyme "
            + enzname
            + ' digests <span class="snppos">'
            + alleletype
            + "</span> allele by "
            + '<span class="enzpos">'
            + enzseq
            + '</span>'
            + ' with mismatches '
            + fmisses(snppos, mismatches)
            + ' .';
    }
    else
    {
        document.getElementById("selectedre").innerHTML = "Enzyme "
            + enzname
            + ' digests <span class="snppos">'
            + alleletype
            + "</span> allele by "
            + '<span class="enzpos">'
            + enzseq
            + '</span>'
            + ' without mismatches.';
    }

//    document.getElementById("selectedre").innerHTML =
//        "Selected RE: {0} {1}".format(enzname, enzseq);
};

function fmisses(snppos, mismatches)
{
    var outstr = '';
    if(mismatches) {
            var mis_posns = mismatches.match(/\d+/g);
            function minussnppos(pos)
            {
                if(pos < (snppos-1)){ return pos - snppos + 1;}
                else { return pos - snppos + 2;}
            };
            function sortpos(a, b) {
                if(a < b) { return -1;}
                if(a > b) { return 1;}
                return 0;
            }
            mis_posns = mis_posns.map(minussnppos);
            mis_posns = mis_posns.sort(sortpos);
            var mis_nucs = mismatches.match(/[a-z]/g);
            outstr += mis_posns.length + ' @ ' + mis_posns.join(', ');
    }
    return outstr;
};