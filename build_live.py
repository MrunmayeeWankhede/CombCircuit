"""Generate the offline live workbench from the verified source anatomy."""
import json
import numpy as np
from combcircuit import ROOT, build, load, verify_hashes

def payload():
    verify_hashes()
    m=build();cells,links,_,_=load(); ids=sorted(cells)
    ci={s:i for i,s in enumerate(ids)}
    edges=[]; aliases={}; canonical=[None]*len(m['C'])
    for sid,g in m['geometry'].items():
        off=int(m['cell_indices'][sid][0])
        for node,local in g['node_index'].items():
            aliases[str(node)]=off+local
            if canonical[off+local] is None:canonical[off+local]=node
        for i,j,length in g['edges']:
            edges.append([off+i,off+j,float(np.pi*.25**2*1e5/(150*length))])
    central=int(m['cell_indices'][4019221][0])
    cut=[central+i for i in m['cut']['local_compartments']]
    cut_edge=next(i for i,e in enumerate(edges) if set(e[:2])==set(cut))
    return dict(version='0.2',baseC=m['C'].tolist(),xyz=[[float(v) for v in p] if np.isfinite(p).all() else None for p in m['xyz']],
        owner=[ci[int(s)] for s in m['owners']],cells=[dict(cells[s],indices=m['cell_indices'][s].tolist()) for s in ids],
        edges=edges,pre=m['pre'].tolist(),post=m['post'].tolist(),feedback=[cells[r['pre']]['type']=='bridge' and cells[r['post']]['type']=='ANN' for r in links],
        aliases=aliases,nodeIds=canonical,stimulus=int(m['stimulus_index']),defaultCut=cut_edge,
        source='Jokura et al. / JekelyLab; see data/provenance.json and docs/DATA_AUDIT.md')

def main():
    d=payload();folder=ROOT/'live'
    (folder/'anatomy.json').write_text(json.dumps(d,separators=(',',':')))
    html=(folder/'template.html').read_text()
    for key,value in [('ANATOMY',json.dumps(d,separators=(',',':'))),('KERNEL',(folder/'solver.js').read_text()),('WORKER',(folder/'worker.js').read_text()),('APP',(folder/'app.js').read_text())]:
        html=html.replace('__'+key+'__',value.replace('</script','<\\/script') if key in ['ANATOMY','KERNEL','WORKER'] else value)
    (ROOT/'LIVE.html').write_text(html)
    print('Built LIVE.html:',len(d['baseC']),'compartments;',len(d['pre']),'contacts')
if __name__=='__main__':main()
