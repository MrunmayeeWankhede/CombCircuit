import unittest
from dataclasses import replace
import numpy as np
from combcircuit import load, verify_hashes, build, simulate, Parameters, arbor_geometry

class ScientificChecks(unittest.TestCase):
    def test_provenance(self):
        self.assertEqual(verify_hashes(),20)
    def test_directed_contacts_and_counts(self):
        c,l,a,d=load()
        self.assertEqual(d['cell_types'],{'balancer':129,'ANN':3,'bridge':26})
        self.assertEqual(len(l),280)
        self.assertEqual(d['balancer_outgoing_contacts'],0)
        self.assertEqual(d['post_only_records'],52)
        self.assertTrue(any(r['pre']==r['post']==4019221 for r in l))
    def test_zero_length_merge(self):
        self.assertEqual(sum(arbor_geometry(a)['merged_zero_length'] for a in load()[2].values()),3)
    def test_capacitance_and_charge(self):
        a,b=build('point'),build('spatial')
        for sid in a['cells']:
            self.assertAlmostEqual(a['C'][a['cell_indices'][sid]].sum(),b['C'][b['cell_indices'][sid]].sum(),places=10)
        self.assertTrue(np.allclose(np.asarray(b['L'].sum(axis=1)).ravel(),0,atol=1e-8))
    def test_zero_input(self):
        r=simulate(build('spatial',replace(Parameters(),pulse_pA=0)))
        self.assertEqual(float(np.max(np.abs(r['traces']))),0)
    def test_passive_cell_average(self):
        p=replace(Parameters(),synaptic_current_pA=0)
        a,b=simulate(build('point',p)),simulate(build('spatial',p))
        np.testing.assert_allclose(a['traces'][:,:5],b['traces'][:,:5],atol=1e-9)
    def test_high_conduction_limit(self):
        p=replace(Parameters(),axial_scale=10000)
        a,b=simulate(build('point',p)),simulate(build('spatial',p))
        np.testing.assert_allclose(a['traces'][:,:5],b['traces'][:,:5],atol=5e-5,rtol=1e-3)
    def test_time_refinement(self):
        a=simulate(build('spatial'));b=simulate(build('spatial',replace(Parameters(),dt_ms=.25)))
        self.assertLess(np.max(np.abs(a['traces'][:,:5]-b['traces'][::2,:5])),.004)
    def test_perturbations(self):
        a,b=build('spatial'),build('spatial',bridge_feedback=False)
        self.assertEqual(len(a['pre'])-len(b['pre']),44)
        c=build('spatial',cut=True)
        self.assertEqual((a['L']-c['L']).nnz,4)
        np.testing.assert_array_equal(a['C'],c['C'])
        with self.assertRaises(ValueError):build('point',cut=True)

if __name__=='__main__':unittest.main()
