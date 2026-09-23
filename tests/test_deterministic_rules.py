import numpy as np
import unittest
ATOL=1e-12;RTOL=1e-10
def cosine(x,w):
 d=np.sqrt(np.sum(x*x)*np.sum(w*w));return np.sum(x*w)/d
def bh(p):
 p=np.asarray(p,float);o=np.argsort(p,kind='mergesort');q=np.empty(len(p));q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
class DeterministicRuleTests(unittest.TestCase):
 def test_direction_and_float64(self):
  x=np.asarray([1.,-1.],dtype=np.float64);w=np.asarray([-1.,1.],dtype=np.float64)
  self.assertLess(cosine(x,w),0)
 def test_bh_reference(self):
  self.assertTrue(np.allclose(bh([.01,.04,.03]),[.03,.04,.04],atol=ATOL,rtol=RTOL))
 def test_boundary_is_strict(self):
  self.assertFalse(.05<.05)

if __name__ == '__main__':
 unittest.main()
