from engine171b import *
from numba import set_num_threads
from concurrent.futures import ThreadPoolExecutor
set_num_threads(8)

def continuous_scores(X,W):
    out=np.empty((len(X),W.shape[1]),dtype=np.float64)
    wn=np.sqrt(np.sum(W*W,axis=0,dtype=np.float64))
    def chunk(lo):
        z=np.asarray(X[lo:lo+256]);zn=np.sqrt(np.sum(z*z,axis=1,dtype=np.float64))
        out[lo:lo+len(z)]=(z@W)/(zn[:,None]*wn[None,:])
    with ThreadPoolExecutor(max_workers=8) as pool:list(pool.map(chunk,range(0,len(X),256)))
    return out

@njit(parallel=True,cache=False,fastmath=False)
def sparse_cosines(z,ids,weights):
    n,k=ids.shape;out=np.empty((z.shape[0],n),dtype=np.float64)
    for i in prange(z.shape[0]):
        for j in range(n):
            numerator=0.;zn=0.;wn=0.
            for g in range(k):
                v=z[i,ids[j,g]];w=weights[j,g]
                numerator+=v*w;zn+=v*v;wn+=w*w
            out[i,j]=numerator/(math.sqrt(zn)*math.sqrt(wn)) if zn>0 and wn>0 else np.nan
    return out

def sparse_scores(X,W):
    ids=np.stack([np.flatnonzero(W[:,j]) for j in range(W.shape[1])])
    weights=np.stack([W[ids[j],j] for j in range(W.shape[1])])
    out=np.empty((len(X),W.shape[1]),dtype=np.float64)
    for lo in range(0,len(X),256):out[lo:lo+256]=sparse_cosines(X[lo:lo+256],ids,weights)
    return out

@njit(parallel=True,cache=False,fastmath=False)
def parallel_median(values,index,ptr):
    out=np.full((len(ptr)-1,values.shape[1]),np.nan)
    for g in prange(len(ptr)-1):
        n=ptr[g+1]-ptr[g]
        if n==1:
            out[g,:]=values[index[ptr[g]],:]
            continue
        for j in range(values.shape[1]):
            v=np.empty(n,dtype=np.float64);k=0
            for a in range(ptr[g],ptr[g+1]):
                x=values[index[a],j]
                if np.isfinite(x):v[k]=x;k+=1
            if k:
                v=np.sort(v[:k]);m=k//2
                out[g,j]=v[m] if k%2 else (v[m-1]+v[m])/2
    return out

def parallel_aggregate(values,h):
    c=parallel_median(values,*h['profile_to_context']);b=parallel_median(c,*h['context_to_base']);p=parallel_median(b,*h['base_to_compound'])
    return c,b,p
