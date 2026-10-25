import numpy as np, time
from epileptor import simulate
SEIZ = 0.0  # x1 above this = ictal-like
def seizing(x1): return (x1 > SEIZ).any(axis=0)  # (trials,N)
t=time.time()
z1 = np.zeros((1,1))
a = simulate(z1, -1.6, T=6000)
b = simulate(z1, -2.4, T=6000)
print("1 isolated x0=-1.6 max x1:", a.max(), "seizes:", seizing(a))
print("2 isolated x0=-2.4 max x1:", b.max(), "seizes:", seizing(b))
# toy 6-region network: region 0 = EZ, strong tie to 1, weak to 2, none to 3-5
C = np.zeros((6,6))
def link(i,j,w): C[i,j]=C[j,i]=w
link(0,1,1.0); link(0,2,0.1); link(3,4,0.5); link(4,5,0.5); link(2,3,0.1)
x0 = np.array([-1.6]+[-2.4]*5)
for K in [0.0, 0.5, 1.0, 2.0, 4.0]:
    o = simulate(C, x0, K=K, T=6000, n_trials=4)
    print("K=",K,"fraction of trials each region seizes:", seizing(o).mean(0).round(2))
print("time", round(time.time()-t,1),"s")
