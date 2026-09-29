"""Full-state agreement between the browser kernel and independent SciPy solves.
Requires Node.js for the JavaScript kernel; the app itself needs only a browser.
"""
import unittest, tempfile, json, subprocess, shutil
from pathlib import Path
import numpy as np
from scipy.sparse import diags, coo_matrix
from scipy.sparse.linalg import factorized
from dataclasses import replace
from combcircuit import ROOT, Parameters, build
from build_live import payload
MAP={'radius':'radius_um','resistivity':'axial_resistivity_ohm_cm','cm':'capacitance_uF_cm2','tau':'leak_tau_ms','cellC':'point_capacitance_pF','gain':'synaptic_current_pA','scale':'synaptic_scale_mV','sign':'synaptic_sign','axial':'axial_scale'}

class LiveAgreement(unittest.TestCase):
 @unittest.skipUnless(shutil.which('node'),'Node.js required to compare browser kernel')
 def test_full_state_and_live_interventions(self):
  d=payload();bridge=next(i for i,c in enumerate(d['cells']) if c['type']=='bridge')
  cases=[{'name':'default','params':{}},{'name':'feedback_off','params':{'feedback':False}},
   {'name':'balanced_cut','params':{'cut':d['defaultCut']}},
   {'name':'arbitrary_cut','params':{'cut':0}},
   {'name':'negative_input','params':{},'amplitude':-1},
   {'name':'negative_thin','params':{'radius':.1,'gain':.1,'sign':-1}},
   {'name':'zero_input','params':{},'amplitude':0},
   {'name':'individual_bridge_input','params':{},'stimulus':d['cells'][bridge]['indices'][0]},
   {'name':'live_change_at_40ms','params':{},'change':{'radius':.5,'gain':.1,'feedback':False,'cut':d['defaultCut']}}]
  for c in cases:
   for k,v in dict(stimulus=d['stimulus'],amplitude=1,duration=20).items():c.setdefault(k,v)
  with tempfile.TemporaryDirectory() as tmp:
   paths=[Path(tmp)/n for n in ['anatomy.json','cases.json','actual.json']]
   paths[0].write_text(json.dumps(d));paths[1].write_text(json.dumps(cases))
   subprocess.run(['node',str(ROOT/'tests/check_live.cjs'),*map(str,paths)],check=True)
   actual=json.loads(paths[2].read_text())
  max_error=0
  for case,out in zip(cases,actual):
   for kind,key in [('spatial','spatial'),('point','point')]:
    settings=case['params'].copy();m=self.model(kind,settings,d);v=np.zeros(len(m['C']));solve=self.solver(m)
    stimulus=case['stimulus'] if kind=='spatial' else d['owner'][case['stimulus']]
    for step in range(320):
     if case.get('change') and step==80:
      settings.update(case['change']);m=self.model(kind,settings,d);solve=self.solver(m)
     p=m['params'];rhs=m['C']/.5*v+np.bincount(m['post'],weights=p.synaptic_sign*p.synaptic_current_pA*np.tanh(v[m['pre']]/p.synaptic_scale_mV),minlength=len(v))
     if 20<=step*.5<20+case['duration']:rhs[stimulus]+=case['amplitude']
     v=solve(rhs)
    err=float(np.max(np.abs(np.asarray(out[key])-v)));max_error=max(max_error,err)
    with self.subTest(case=case['name'],model=kind):self.assertLess(err,2e-8)
  print('Live JS vs SciPy maximum full-state error:',max_error,'mV across',len(cases),'scenarios')
 def model(self,kind,s,d):
  p=replace(Parameters(),**{MAP[k]:v for k,v in s.items() if k in MAP})
  m=build(kind,p,bridge_feedback=s.get('feedback',True))
  if kind=='spatial' and s.get('cut',-1)>=0:
   a,b,g=d['edges'][s['cut']];g*= (p.radius_um/.25)**2*(150/p.axial_resistivity_ohm_cm)*p.axial_scale
   m['L']=m['L']-coo_matrix(([g,g,-g,-g],([a,b,a,b],[a,b,b,a])),shape=m['L'].shape).tocsc()
  return m
 def solver(self,m):return factorized(diags(m['C']/.5+m['C']/m['params'].leak_tau_ms)+m['L'])
if __name__=='__main__':unittest.main()
