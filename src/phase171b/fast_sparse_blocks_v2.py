from engine171b import *
from numba import set_num_threads
set_num_threads(8)

@njit(parallel=True,cache=False,fastmath=False)
def sparse_blocks_256(z,ids,weights):
    n,k=ids.shape;out=np.empty((z.shape[0],n),dtype=np.float64)
    for block in prange((z.shape[0]+255)//256):
        for i in range(block*256,min((block+1)*256,z.shape[0])):
            for j in range(n):
                numerator=0.;zn=0.;wn=0.
                for g in range(k):
                    v=z[i,ids[j,g]];w=weights[j,g]
                    numerator+=v*w;zn+=v*v;wn+=w*w
                out[i,j]=numerator/(math.sqrt(zn)*math.sqrt(wn)) if zn>0 and wn>0 else np.nan
    return out

def sparse_scores_v2(X,W):
    ids=np.stack([np.flatnonzero(W[:,j]) for j in range(W.shape[1])])
    weights=np.stack([W[ids[j],j] for j in range(W.shape[1])])
    return sparse_blocks_256(X,ids,weights)
