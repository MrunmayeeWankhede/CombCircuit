"""A reproducible, exploratory model of measured comb-jelly circuit anatomy.

Run: python combcircuit.py
No biological parameters are inferred by this prototype. See docs/MODEL.md.
"""
from pathlib import Path
from dataclasses import dataclass, asdict, replace
from collections import Counter
import csv
import json
import hashlib
import argparse
import platform
import numpy as np
import scipy
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import factorized

ROOT = Path(__file__).resolve().parent
UPSTREAM = ROOT / 'data/upstream'
ANN_NAMES = {4019221: 'ANN Q1-4', 4018492: 'ANN Q1Q2', 4018688: 'ANN Q3Q4'}

def rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))

def verify_hashes():
    manifest = json.loads((ROOT/'data/provenance.json').read_text())
    for f in manifest['files']:
        if hashlib.sha256((UPSTREAM/f['path']).read_bytes()).hexdigest() != f['sha256']:
            raise ValueError('Checksum mismatch: '+f['path'])
    return len(manifest['files'])

def load():
    master = rows(UPSTREAM/'analysis/data/stats_master.csv')
    cells = {int(r['skid']): {'id':int(r['skid']), 'type':'ANN' if r['celltype']=='SSN' else r['celltype'],
                              'name':ANN_NAMES.get(int(r['skid']),r['celltype']+' '+r['skid'])}
             for r in master if r['celltype'] in {'SSN','bridge','balancer'}}
    contacts = rows(UPSTREAM/'analysis/data/stats_synapse.csv')
    pre, post = {}, {}
    for raw in contacts:
        r = {k:int(raw[k]) for k in ['treenode_id','connector_id','prepost','skid']}
        cid = r['connector_id']
        if r['prepost']==0:
            if cid in pre:
                raise ValueError('Multiple presynaptic endpoints')
            pre[cid]=r
        elif r['prepost']==1:
            post.setdefault(cid,[]).append(r)
        else:
            raise ValueError('Unexpected contact relation')
    links=[]
    all_matched=0
    for cid,p in pre.items():
        for q in post.get(cid,[]):
            all_matched+=1
            if p['skid'] in cells and q['skid'] in cells:
                links.append({'connector_id':cid,'pre':p['skid'],'post':q['skid'],
                              'pre_node':p['treenode_id'],'post_node':q['treenode_id']})
    arbors={sid:json.loads((UPSTREAM/f'catmaid/{sid}.json').read_text()) for sid in ANN_NAMES}
    for sid,a in arbors.items():
        ids={r[0] for r in a[0]}
        required={r[key] for r in links for key,side in [('pre_node','pre'),('post_node','post')] if r[side]==sid}
        if required-ids:
            raise ValueError(f'Source synaptic treenodes missing from arbor {sid}')
    audit={'scope':'ANN/bridge/balancer induced subgraph', 'cells':len(cells),
           'cell_types':dict(Counter(c['type'] for c in cells.values())),
           'raw_contact_rows':len(contacts),'unique_connectors':len(set(pre)|set(post)),
           'presynaptic_connectors':len(pre),'matched_links_all_types':all_matched,
           'post_only_connectors':sorted(set(post)-set(pre)),
           'post_only_records':sum(len(post[c]) for c in set(post)-set(pre)),
           'pre_only_connectors':sorted(set(pre)-set(post)),
           'focus_contacts':len(links),'focus_directed_pairs':len({(r['pre'],r['post']) for r in links}),
           'balancer_outgoing_contacts':sum(cells[r['pre']]['type']=='balancer' for r in links),
           'figure3_group_contact_sum':sum(int(r['Freq']) for r in rows(UPSTREAM/'manuscript/source_data/Figure3_source_data2.csv')),
           'raw_arbor_nodes':{str(s):len(a[0]) for s,a in arbors.items()},
           'notes':['SSN source label mapped to ANN; SNN is a separate label.',
                    'Orphan records excluded without inventing missing partners.',
                    'Live CATMAID snapshot and pinned GitHub contacts are distinct versions; contact IDs checked.',
                    'Rooted skeleton trees do not establish all internal syncytial loops.',
                    'Raw CATMAID connectors retained but not mixed into the pinned contact table.']}
    if audit['focus_contacts'] != audit['figure3_group_contact_sum']:
        raise ValueError('Figure 3 contact total mismatch')
    return cells, links, arbors, audit

@dataclass(frozen=True)
class Parameters:
    radius_um: float = .25
    axial_resistivity_ohm_cm: float = 150.
    capacitance_uF_cm2: float = 1.
    leak_tau_ms: float = 20.
    point_capacitance_pF: float = .5
    synaptic_current_pA: float = .05
    synaptic_scale_mV: float = 5.
    synaptic_sign: int = 1
    axial_scale: float = 1.
    pulse_pA: float = 1.
    pulse_start_ms: float = 20.
    pulse_end_ms: float = 40.
    duration_ms: float = 160.
    dt_ms: float = .5

    def validate(self):
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Parameters must be finite')
        for key in ['radius_um','axial_resistivity_ohm_cm','capacitance_uF_cm2','leak_tau_ms','point_capacitance_pF','synaptic_scale_mV','dt_ms','duration_ms']:
            if getattr(self,key)<=0:
                raise ValueError(key+' must be positive')
        if self.synaptic_sign not in [-1,1] or self.axial_scale<0 or self.synaptic_current_pA<0:
            raise ValueError('Invalid sign or scale')
        if not 0<=self.pulse_start_ms<self.pulse_end_ms<=self.duration_ms:
            raise ValueError('Pulse must lie inside simulation')
        for t in [self.pulse_start_ms,self.pulse_end_ms,self.duration_ms]:
            if abs(t/self.dt_ms-round(t/self.dt_ms))>1e-7:
                raise ValueError('Pulse and duration must align with timestep')

def arbor_geometry(arbor):
    raw=arbor[0]
    ids=[n[0] for n in raw]
    idx={n:i for i,n in enumerate(ids)}
    xyz=np.array([n[3:6] for n in raw],dtype=float)/1000
    parent=np.arange(len(ids))
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]]
            i=parent[i]
        return i
    edges=[]
    for n in raw:
        if n[1] is None:
            continue
        if n[1] not in idx:
            raise ValueError('Missing arbor parent')
        i,j=idx[n[0]],idx[n[1]]
        length=np.linalg.norm(xyz[i]-xyz[j])
        if length<1e-9:
            parent[root(i)]=root(j)
        else:
            edges.append((i,j,length))
    roots=[root(i) for i in range(len(ids))]
    unique=sorted(set(roots));reduced={r:i for i,r in enumerate(unique)}
    oldnew=np.array([reduced[r] for r in roots])
    edges=[(int(oldnew[i]),int(oldnew[j]),length) for i,j,length in edges if oldnew[i]!=oldnew[j]]
    lengths=np.zeros(len(unique))
    for i,j,length in edges:
        lengths[i]+=length/2;lengths[j]+=length/2
    if (lengths<=0).any():
        raise ValueError('Isolated compartment')
    return {'xyz':xyz[unique],'edges':edges,'lengths':lengths,
            'node_index':{nid:int(oldnew[i]) for i,nid in enumerate(ids)},
            'merged_zero_length':len(ids)-len(unique)}

def build(kind='spatial',params=None,bridge_feedback=True,cut=False):
    p=params or Parameters();p.validate()
    if kind not in {'spatial','point'}:
        raise ValueError('Unknown model')
    if cut and kind=='point':
        raise ValueError('A point model cannot represent a branch cut')
    cells,links,arbors,audit=load()
    geometry={sid:arbor_geometry(a) for sid,a in arbors.items()}
    C,owners,xyz,contact_ix,cell_indices=[],[],[],{},{}
    cable_edges,offset=[],{}
    for sid in sorted(cells):
        start=len(C);offset[sid]=start
        if sid in geometry:
            g=geometry[sid]
            cap=.01*p.capacitance_uF_cm2*2*np.pi*p.radius_um*g['lengths']
            if kind=='spatial':
                C.extend(cap);xyz.extend(g['xyz']);owners.extend([sid]*len(cap))
                for n,i in g['node_index'].items():contact_ix[(sid,n)]=start+i
                for i,j,length in g['edges']:
                    cond=p.axial_scale*np.pi*p.radius_um**2*1e5/(p.axial_resistivity_ohm_cm*length)
                    cable_edges.append((start+i,start+j,cond,sid))
            else:
                C.append(cap.sum());xyz.append(np.average(g['xyz'],axis=0,weights=cap));owners.append(sid)
                for n in g['node_index']:contact_ix[(sid,n)]=start
        else:
            C.append(p.point_capacitance_pF);xyz.append([np.nan]*3);owners.append(sid)
        cell_indices[sid]=np.arange(start,len(C))
    C,owners,xyz=np.array(C),np.array(owners),np.array(xyz)
    central=4019221;g=geometry[central]
    contact_nodes=sorted({r[key] for r in links for key,side in [('pre_node','pre'),('post_node','post')] if r[side]==central})
    stimulus_node=min(contact_nodes,key=lambda n:(g['xyz'][g['node_index'][n],0],n))
    stimulus_index=contact_ix[(central,stimulus_node)]
    # Pick a reproducible, balanced-length tree cut. This is not a biological lesion claim.
    adj=[[] for _ in g['lengths']]
    for k,(i,j,l) in enumerate(g['edges']):
        adj[i].append((j,k));adj[j].append((i,k))
    traversal,parents,edge_parent=[0],{0:-1},{}
    for i in traversal:
        for j,k in adj[i]:
            if j not in parents:
                parents[j]=i;edge_parent[j]=k;traversal.append(j)
    if len(traversal)!=len(adj) or len(g['edges'])!=len(adj)-1:
        raise ValueError('Expected central export to be a connected tree')
    mass=g['lengths'].copy()
    for i in reversed(traversal[1:]):mass[parents[i]]+=mass[i]
    child=min(traversal[1:],key=lambda i:abs(mass[i]/mass[0]-.5))
    cut_local=g['edges'][edge_parent[child]][:2]
    cut_global={offset[central]+i for i in cut_local}
    cut_info={'cell':central,'local_compartments':list(cut_local),'partition_fraction':float(mass[child]/mass[0]),
              'selection':'edge splitting traced cable length most evenly'}
    rr,cc,vv=[],[],[]
    for i,j,cond,sid in cable_edges:
        if cut and sid==central and {i,j}==cut_global:continue
        rr.extend([i,j,i,j]);cc.extend([i,j,j,i]);vv.extend([cond,cond,-cond,-cond])
    L=coo_matrix((vv,(rr,cc)),shape=(len(C),len(C))).tocsc()
    pre,post=[],[]
    for r in links:
        if not bridge_feedback and cells[r['pre']]['type']=='bridge' and cells[r['post']]['type']=='ANN':continue
        pre.append(contact_ix[(r['pre'],r['pre_node'])] if r['pre'] in geometry else int(cell_indices[r['pre']][0]))
        post.append(contact_ix[(r['post'],r['post_node'])] if r['post'] in geometry else int(cell_indices[r['post']][0]))
    return dict(kind=kind,params=p,cells=cells,C=C,owners=owners,xyz=xyz,cell_indices=cell_indices,L=L,
                pre=np.array(pre),post=np.array(post),stimulus_index=stimulus_index,stimulus_node=stimulus_node,
                cut=cut_info,geometry=geometry,audit=audit,bridge_feedback=bridge_feedback,cut_enabled=cut)

def simulate(model,snapshot_count=33):
    p,C=model['params'],model['C']
    steps=round(p.duration_ms/p.dt_ms);t=np.arange(steps+1)*p.dt_ms
    solve=factorized(diags(C/p.dt_ms+C/p.leak_tau_ms)+model['L'])
    v=np.zeros(len(C))
    groups=[(name,model['cell_indices'][sid]) for sid,name in ANN_NAMES.items()]
    types=np.array([model['cells'][int(s)]['type'] for s in model['owners']])
    groups += [(typ,np.where(types==typ)[0]) for typ in ['bridge','balancer']]
    labels=[n for n,i in groups]+['stimulated compartment']
    traces=np.zeros((len(t),len(labels)))
    snap_steps=set(np.linspace(0,steps,snapshot_count,dtype=int))
    ix=np.where(np.isfinite(model['xyz'][:,0]))[0]
    snapshots,snapshot_t=[v[ix].copy()],[0.]
    for step in range(steps):
        syn=np.bincount(model['post'],weights=p.synaptic_sign*p.synaptic_current_pA*np.tanh(v[model['pre']]/p.synaptic_scale_mV),minlength=len(v))
        rhs=C/p.dt_ms*v+syn
        if p.pulse_start_ms<=t[step]<p.pulse_end_ms:rhs[model['stimulus_index']]+=p.pulse_pA
        v=solve(rhs)
        if not np.isfinite(v).all():raise FloatingPointError('Nonfinite voltage')
        for k,(_,indices) in enumerate(groups):traces[step+1,k]=np.average(v[indices],weights=C[indices])
        traces[step+1,-1]=v[model['stimulus_index']]
        if step+1 in snap_steps:snapshots.append(v[ix].copy());snapshot_t.append(float(t[step+1]))
    metrics={}
    for k,label in enumerate(labels):
        peak=int(np.argmax(np.abs(traces[:,k])))
        integral=np.trapezoid(traces[:,k],t) if hasattr(np,'trapezoid') else np.trapz(traces[:,k],t)
        metrics[label]={'peak_abs_mV':float(abs(traces[peak,k])),'peak_time_ms':float(t[peak]),'integral_mV_ms':float(integral)}
    return dict(t=t,labels=labels,traces=traces,snapshots=np.array(snapshots),snapshot_t=snapshot_t,snapshot_indices=ix,metrics=metrics)

def export_csv(path,fields,records):
    with open(path,'w',newline='') as f:
        w=csv.writer(f);w.writerow(fields);w.writerows(records)

def make_figure(out,models,results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,3,figsize=(14,4.5))
    for (sid,g),color in zip(models['spatial']['geometry'].items(),['#167d99','#d06682','#dc9a28']):
        segments=[[g['xyz'][i,[0,1]],g['xyz'][j,[0,1]]] for i,j,l in g['edges']]
        axs[0].add_collection(LineCollection(segments,colors=color,linewidths=.45,alpha=.8,label=ANN_NAMES[sid]))
    axs[0].autoscale();axs[0].set_aspect('equal');axs[0].set_xlabel('x (µm)');axs[0].set_ylabel('y (µm)')
    axs[0].set_title('Measured ANN traces — xy projection');axs[0].legend(fontsize=8)
    for name,color,style in [('point','#617084','--'),('spatial','#167d99','-'),('no_bridge_feedback','#d06682','-'),('branch_cut','#dc9a28','-')]:
        r=results[name]
        axs[1].plot(r['t'],r['traces'][:,0],style,color=color,label=name.replace('_',' '))
        axs[2].plot(r['t'],r['traces'][:,4],style,color=color)
    for ax in axs[1:]:
        ax.axvspan(20,40,color='#c6dce0',alpha=.35);ax.set_xlabel('Time (ms)');ax.set_ylabel('Mean deviation from rest (mV)')
    axs[1].set_title('Central ANN — capacitance-weighted mean');axs[1].legend(fontsize=8)
    axs[2].set_title('Balancer cells — modeled electrical drive')
    fig.suptitle('CombCircuit 0.1 • hypothetical physiology, real circuit anatomy',fontsize=14,y=1.01)
    fig.text(.5,-.02,'Illustrative pulse response. No measured neural voltage, ciliary beat prediction or validation claim.',ha='center',fontsize=10)
    fig.tight_layout();fig.savefig(out/'comparison.png',dpi=170,bbox_inches='tight');plt.close(fig)

def make_viewer(path,models,results,audit,summary):
    m=models['spatial'];ix=results['spatial']['snapshot_indices'];inverse={int(v):i for i,v in enumerate(ix)}
    edges=[]
    for sid,g in m['geometry'].items():
        start=int(m['cell_indices'][sid][0]);edges += [[inverse[start+i],inverse[start+j]] for i,j,l in g['edges']]
    owners=m['owners'][ix]
    payload={'xyz':np.round(m['xyz'][ix],3).tolist(),'owners':owners.tolist(),'edges':edges,'names':ANN_NAMES,
             'audit':audit,'summary':summary,'stimulus':inverse[m['stimulus_index']],
             'cut_edge':[inverse[int(m['cell_indices'][4019221][0])+i] for i in m['cut']['local_compartments']],
             'matrix':rows(UPSTREAM/'manuscript/source_data/Figure3_source_data2.csv'),'cases':{}}
    for name,r in results.items():
        snapshots=r['snapshots']
        if name=='point':
            lookup={int(s):i for i,s in enumerate(models[name]['owners'][r['snapshot_indices']])}
            snapshots=snapshots[:,[lookup[int(s)] for s in owners]]
        payload['cases'][name]={'t':r['t'].tolist(),'labels':r['labels'],'traces':np.round(r['traces'],6).tolist(),
                                'snapshot_t':r['snapshot_t'],'voltage':np.round(snapshots,4).tolist()}
    data=json.dumps(payload,separators=(',',':'),allow_nan=False).replace('</','<\\/')
    path.write_text((ROOT/'viewer_template.html').read_text().replace('__PAYLOAD__',data))

def run(out=None,quick=False):
    out=Path(out or ROOT/'results');out.mkdir(parents=True,exist_ok=True)
    verify_hashes();cells,links,arbors,audit=load()
    (out/'audit.json').write_text(json.dumps(audit,indent=2))
    params=Parameters();models,results={},{}
    for name,kind,feedback,cut in [('point','point',True,False),('spatial','spatial',True,False),('no_bridge_feedback','spatial',False,False),('branch_cut','spatial',True,True)]:
        m=build(kind,params,bridge_feedback=feedback,cut=cut);r=simulate(m)
        models[name],results[name]=m,r
        export_csv(out/f'{name}_traces.csv',['time_ms']+r['labels'],zip(r['t'],*r['traces'].T))
        print(f'{name}: {len(m["C"]):,} compartments, {len(m["pre"])} contacts')
    export_csv(out/'contacts.csv',list(links[0]),[list(r.values()) for r in links])
    sweep=[]
    if not quick:
        for radius in [.1,.25,.5]:
            for sign in [-1,1]:
                for gain in [0,.05,.1]:
                    p=replace(params,radius_um=radius,synaptic_sign=sign,synaptic_current_pA=gain)
                    a,b=simulate(build('point',p)),simulate(build('spatial',p))
                    sweep.append([radius,sign,gain,float(np.sqrt(np.mean((a['traces'][:,:5]-b['traces'][:,:5])**2))),a['metrics']['balancer']['peak_abs_mV'],b['metrics']['balancer']['peak_abs_mV']])
        export_csv(out/'sensitivity.csv',['radius_um','all_synapses_sign','current_per_contact_pA','mean_trace_rmse_mV','point_balancer_peak_mV','spatial_balancer_peak_mV'],sweep)
    summary={'status':'Exploratory synthetic pulse response; no physiological fit or novelty claim',
             'parameters':asdict(params),'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
             'stimulus_treenode':models['spatial']['stimulus_node'],'branch_cut':models['spatial']['cut'],
             'metrics':{n:r['metrics'] for n,r in results.items()},'sensitivity_scenarios':len(sweep)}
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    make_figure(out,models,results);make_viewer(out.parent/'EXPLORE.html',models,results,audit,summary)
    print('Open',out.parent/'EXPLORE.html','in a browser; no server or internet needed.')
    return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',default=str(ROOT/'results'))
    parser.add_argument('--quick',action='store_true',help='Skip sensitivity sweep')
    a=parser.parse_args();run(a.out,a.quick)
