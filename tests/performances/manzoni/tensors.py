import timeit

setup = """
import numpy as np
x = np.random.rand(int(1e7),4)
T = np.random.rand(4,4,4)
"""

s = """
Tx = np.tensordot(T,x,axes=((1),(1)))
out = np.einsum('ikl,lk->li',Tx,x)
"""
t = timeit.timeit(s, setup=setup, number=1)
print("====================================================")
print(f"Solution 1:")
print(f"{s}")
print(f"===> {1000*t} ms")

s = """
np.einsum('ijl,lj->li', T@x.T, x)
"""
t = timeit.timeit(s, setup=setup, number=1)
print("====================================================")
print(f"Solution 2:")
print(f"{s}")
print(f"===> {1000*t} ms")

# https://stackoverflow.com/a/51048168/420892
s = """
W = np.matmul(T,x.T)
ZT = np.sum(W*x.T[np.newaxis,:,:], axis=1)
Z = ZT.T
"""
t = timeit.timeit(s, setup=setup, number=1)
print("====================================================")
print(f"Solution 3:")
print(f"{s}")
print(f"===> {1000*t} ms")

# https://stackoverflow.com/a/51047506/420892
s = """
np.einsum('ijk,lj,lk->li', T, x, x)
"""
t = timeit.timeit(s, setup=setup, number=1)
print("====================================================")
print(f"Solution 3:")
print(f"{s}")
print(f"===> {1000*t} ms")

s = """
((x[:, None, None, :]@T).squeeze()@x[..., None]).squeeze()
"""
t = timeit.timeit(s, setup=setup, number=1)
print("====================================================")
print(f"Solution 3:")
print(f"{s}")
print(f"===> {1000*t} ms")

s = """
np.einsum('ijl,lj->li', T@x.T, x)
"""
t = timeit.timeit(s, setup=setup, number=1)
print("====================================================")
print(f"Solution 3:")
print(f"{s}")
print(f"===> {1000*t} ms")

# https://stackoverflow.com/a/51048025/420892
s = """
Tx = np.tensordot(T,x,axes=((1),(1)))
out = np.einsum('ikl,lk->li',Tx,x)
"""
t = timeit.timeit(s, setup=setup, number=1)
print("====================================================")
print(f"Solution 3:")
print(f"{s}")
print(f"===> {1000*t} ms")


setup = """
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
"""

s = """
res=tensor_mult(T,x)
"""

t = timeit.timeit(s, setup=setup, number=1)
print("====================================================")
print(f"Solution 3:")
print(f"{s}")
print(f"===> {1000*t} ms")
