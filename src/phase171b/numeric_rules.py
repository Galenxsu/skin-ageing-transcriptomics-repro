import math
ATOL=1e-12
RTOL=1e-10
def state(x):
    if x is None or isinstance(x,str) and x=='NA':return 'NA'
    if math.isnan(float(x)):return 'NaN'
    if math.isinf(float(x)):return 'Inf' if x>0 else '-Inf'
    return 'FINITE'
def direction(x):return 'POSITIVE' if x>0 else 'NEGATIVE' if x<0 else 'ZERO'
def compare(x1,x2,decision=None):
    a,b=state(x1),state(x2)
    r=dict(comparison_value=str(x1),reference_value=str(x2),state1=a,state2=b,absolute_difference='NA',relative_difference='NA',numeric_pass=False,decision_pass=False,status='MISSING_STATE_DISCORDANCE')
    if a!=b:return r
    if a!='FINITE':r.update(numeric_pass=True,decision_pass=True,status='MATCHING_NONFINITE_STATE');return r
    d=abs(x1-x2);rel=d/abs(x2) if x2 else (0. if d==0 else 'Inf')
    ok=d<=ATOL+RTOL*abs(x2);r.update(absolute_difference=d,relative_difference=rel,numeric_pass=ok,decision_pass=True,status='PASS' if ok else 'NUMERICAL_TOLERANCE_FAILURE')
    if decision=='direction' and direction(x1)!=direction(x2):r.update(decision_pass=False,status='DIRECTION_DISCORDANCE')
    if decision=='significance' and (x1<.05)!=(x2<.05):r.update(decision_pass=False,status='DECISION_BOUNDARY_DISCORDANCE')
    return r
def unit_tests():
    cases=[(1.,1.+1e-13,None,'PASS'),(-1e-15,1e-15,'direction','DIRECTION_DISCORDANCE'),(0.,1e-15,'direction','DIRECTION_DISCORDANCE'),(.05-1e-15,.05+1e-15,'significance','DECISION_BOUNDARY_DISCORDANCE'),(None,float('nan'),None,'MISSING_STATE_DISCORDANCE'),(float('inf'),float('-inf'),None,'MISSING_STATE_DISCORDANCE'),(float('nan'),float('nan'),None,'MATCHING_NONFINITE_STATE'),(float('inf'),float('inf'),None,'MATCHING_NONFINITE_STATE'),(0.,0.,'direction','PASS'),(1.,1.01,None,'NUMERICAL_TOLERANCE_FAILURE')]
    rr=[]
    for a,b,k,expected in cases:
        r=compare(a,b,k);assert r['status']==expected;r['expected']=expected;rr.append(r)
    return rr
