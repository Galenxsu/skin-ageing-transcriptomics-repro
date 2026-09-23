from common171b import *
import numpy as np

maximum_zeros=0;total_zeros=0
for lo in range(0,107201,256):
    source=np.load(O/'MATRIX/frozen_107201x8227_float64.npy',mmap_mode='r')
    counts=np.count_nonzero(source[lo:lo+256]==0,axis=1);maximum_zeros=max(maximum_zeros,int(counts.max()));total_zeros+=int(counts.sum());source._mmap.close()
scores=np.load(O/'RESULTS/observed_profile_scores.npy');names=jr(O/'RESULTS/query_order.json')
rows=[dict(query_id=q,nonfinite_observed_profiles=int(np.count_nonzero(~np.isfinite(scores[:,i])))) for i,q in enumerate(names)]
csvwrite(O/'QA/QUERY_NUMERICAL_INTEGRITY.csv',rows)
js(O/'QA/NULL_ZERO_NORM_STRUCTURAL_CHECK.json',dict(utc=now(),total_exact_zero_entries=total_zeros,maximum_zero_entries_per_profile=maximum_zeros,minimum_nonzero_query_members=50,zero_norm_impossible_for_any_frozen_null_query=maximum_zeros<50,explanation='If each profile has fewer than 50 zero elements among 8227, no permitted query with at least 50 distinct nonzero-weight genes can have zero expression norm. This is an integrity check, not gene selection.'))
print('MAX_ZEROS_PER_PROFILE',maximum_zeros,'NULL_ZERO_NORM_IMPOSSIBLE',maximum_zeros<50,flush=True)
