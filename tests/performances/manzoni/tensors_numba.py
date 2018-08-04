import numba as nb
import numpy as np
import time

@nb.njit(fastmath=True,parallel=True)
def tensor_mult(T,x):
  res=np.empty((x.shape[0],T.shape[0]),dtype=T.dtype)
  for l in nb.prange(x.shape[0]):
    for i in range(T.shape[0]):
      sum=0.
      for j in range(T.shape[1]):
        for k in range(T.shape[2]):
          sum+=T[i,j,k]*x[l,j]*x[l,k]
      res[l,i]=sum
  return res

x = np.random.rand(1000000,6)
T = np.random.rand(6,6,6)

#first call has a compilation overhead (about 0.6s)
res=tensor_mult(T,x)

t1=time.time()
for i in range(10):
  #@divakar
  #Tx = np.tensordot(T,x,axes=((1),(1)))
  #out = np.einsum('ikl,lk->li',Tx,x)

  res=tensor_mult(T,x)

print(time.time()-t1)
