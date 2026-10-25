#!/usr/bin/env python3
"""Check cross-file reported results without third-party packages.

This validates committed reference CSV consistency; it is NOT an independent
recomputation of exact trial probabilities or a substitute for running the R code.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / 'reference_results'

def load(rel):
    with (REF / rel).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def approx(x, expected, tolerance=0.00015):
    assert abs(float(x)-expected) <= tolerance, f'Value {x} is not within {tolerance} of {expected}'

nulls = load('R2_adjusted/R2_C4_Point_global_null_percent.csv')
alts = load('R2_adjusted/R2_C4_Original_30_alternatives_percent.csv')
sub = load('R2_adjusted/R2_C4_Subgroup_only_percent.csv')
caps = load('R2_adjusted/R2_C4_component_cap_search.csv')
assert len(nulls) == 20 and len(alts) == 120 and len(sub) == 8
methods = sorted({r['method'] for r in alts})
assert len(methods) == 4, methods
proposed = next(m for m in methods if m == 'Proposed')
original = next(m for m in methods if m == 'Independent_BOP2_unadjusted')
adjusted = next(m for m in methods if m == 'Independent_BOP2_global_error_matched')

def max_rate(method):
    return max(float(r['PRN_any']) for r in nulls if r['method']==method)

approx(max_rate(proposed), 0.09894542782786911,1e-10)
approx(max_rate(original), 0.12633818863934865,1e-10)
approx(max_rate(adjusted), 0.0989178353869672,1e-10)

for rows in (nulls, alts, sub):
    for r in rows:
        approx(float(r['PRN_any'])-float(r['PRN_all'])-float(r['PRN_positive']), 0, 5e-10)

for method, expected in [(proposed, 0.8345),(original,0.7876),(adjusted,0.7733)]:
    arr=[float(r['PRN_any']) for r in alts if r['method']==method]
    assert len(arr)==30
    approx(sum(arr)/30,expected,0.00007)

by_scenario = lambda m: {(r['pi'],r['theta_positive'],r['theta_negative']): float(r['PRN_any']) for r in alts if r['method']==m}
p = by_scenario(proposed)
a = by_scenario(adjusted)
assert len(p)==30 and len(a)==30
wins=[k for k in p if p[k]>a[k]]
assert len(wins)==29
exception = next(k for k in p if p[k] <= a[k])
assert tuple(map(float,exception))==(0.4,0.6,0.2), exception
approx(p[exception],0.829622,0.000005)
approx(a[exception],0.837452,0.000005)

scenario=lambda method: next(r for r in sub if r['method']==method and float(r['theta_positive'])==0.4)
s=scenario(proposed)
assert abs(float(s['theta_all'])-0.19)<1e-10
approx(s['PRN_positive_pct'],28.72,0.01)
approx(s['PRN_all_pct'],4.63,0.01)

feasible=[r for r in caps if r['valid_nominal_10pct'].lower()=='true']
assert feasible
best = max(feasible, key=lambda r: (float(r['max_point_null']),float(r['component_cap'])))
approx(best['component_cap'],0.0874601406339,1e-12)
approx(best['max_point_null'],0.0989178353869672,1e-12)

r1=load('R1_reversed/R1M2_reverse_table_percent.csv')
assert len(r1)==20
print('PASS: point-global-null rates, claim decomposition, 30-scenario means,')
print('      29-of-30 comparison, subgroup-only case, and 1050-candidate tau selection.')
print('NOTE: no R program or exact enumeration was executed by this check.')

# Reviewer 1 benchmark results: the main manuscript states gamma <= 8.
primary_summary = load('R1_benchmarks/primary_gamma8/ALL_summary.csv')
primary_selected = load('R1_benchmarks/primary_gamma8/ALL_selected_designs.csv')
assert len(primary_summary) == len(primary_selected) == 6
by_key = lambda rows: {(float(r['p0']),float(r['p1']),r['Method']):r for r in rows}
prim_sum = by_key(primary_summary)
prim_sel = by_key(primary_selected)
assert prim_sum.keys() == prim_sel.keys()
assert {(k[0],k[1]) for k in prim_sum} == {(0.1,0.3),(0.2,0.4),(0.3,0.5)}
expected_primary = {
    (0.1,0.3,'Globally calibrated'):(0.0923024945859645,0.835106837222923),
    (0.2,0.4,'Globally calibrated'):(0.098945427827869,0.764425277321558),
    (0.3,0.5,'Globally calibrated'):(0.0979694473765788,0.72609168172837),
    (0.1,0.3,'Independent conventional BOP2 components'):(0.100078823847204,0.784845784740067),
    (0.2,0.4,'Independent conventional BOP2 components'):(0.126338188639349,0.70325350043878),
    (0.3,0.5,'Independent conventional BOP2 components'):(0.146137889107127,0.756352494942012),
}
for k,(alpha0,power0) in expected_primary.items():
    approx(prim_sum[k]['max_point_null_PRN_any'],alpha0,1e-11)
    approx(prim_sum[k]['mean_working_alternative_PRN_any'],power0,1e-11)
    approx(prim_sel[k]['max_global_type1'],alpha0,1e-11)
    approx(prim_sel[k]['average_power'],power0,1e-11)
    if k[2]=='Globally calibrated':
        assert float(prim_sel[k]['gamma_positive'])<=8.0

extended_sum = load('R1_benchmarks/exploratory_gamma12/ALL_summary.csv')
extended_sel = load('R1_benchmarks/exploratory_gamma12/ALL_selected_designs.csv')
assert len(extended_sum) == len(extended_sel) == 4
ext_sum = by_key(extended_sum)
ext_sel = by_key(extended_sel)
assert ext_sum.keys()==ext_sel.keys()
assert {(k[0],k[1]) for k in ext_sum} == {(0.1,0.3),(0.3,0.5)}
for k in ext_sum:
    assert k in prim_sum
    # The two searches differ at high numerical precision, but not at
    # the displayed two-decimal percentage precision of Table S7.
    assert abs(float(ext_sum[k]['max_point_null_PRN_any'])-float(prim_sum[k]['max_point_null_PRN_any']))<5e-8
    assert abs(float(ext_sum[k]['mean_working_alternative_PRN_any'])-float(prim_sum[k]['mean_working_alternative_PRN_any']))<5e-8
    if k[2]=='Globally calibrated':
        assert float(ext_sel[k]['gamma_positive'])>8.0

# Reproduction of the published BOP2 Table I boundaries and same-rule
# component recalibration to the present paper's 20/30/40 schedule.
published = load('BOP2_validation/bop2_reference_reproduction.csv')
allcomer_bop2 = load('BOP2_validation/bop2_allcomer_current_schedule.csv')
positive_library = load('BOP2_validation/bop2_positive_entry_specific_library.csv')
positive_boundaries = load('BOP2_validation/bop2_positive_entry_specific_boundaries.csv')
assert len(published)==len(allcomer_bop2)==1
assert published[0]['Looks']=='10,15,20,25,30,35,40'
assert published[0]['Boundaries']=='1,2,4,5,7,9,10'
approx(published[0]['exact_type1'],0.0974284191403399,1e-12)
approx(published[0]['exact_power'],0.881487130225181,1e-12)
assert allcomer_bop2[0]['Looks']=='20,30,40'
assert allcomer_bop2[0]['Boundaries']=='5,7,10'
approx(allcomer_bop2[0]['exact_type1'],0.0949197824332817,1e-12)
approx(allcomer_bop2[0]['exact_power'],0.860683798765972,1e-12)
assert [int(r['first_assessment_n']) for r in positive_library]==list(range(10,31))
assert len(positive_boundaries)>21

print('PASS: supplied R1 benchmark gamma-8/12 snapshots and published BOP2')
print('      Table I stopping-boundary reproduction snapshots are consistent.')
