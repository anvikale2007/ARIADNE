import numpy as np, time, warnings
warnings.filterwarnings("ignore")
from epileptor_fast import simulate
from connectome import load_template
from groups import coarse, GROUP_NAMES
t = time.time()
a = simulate(np.zeros((1,1)), -1.6); b = simulate(np.zeros((1,1)), -2.4)
print("1 isolated -1.6 seizes:", bool((a>0).any()), "| 2 isolated -2.4 seizes:", bool((b>0).any()), "| first call (compile) s:", round(time.time()-t,1))
C = np.zeros((6,6))
def link(i,j,w): C[i,j]=C[j,i]=w
link(0,1,1.0); link(0,2,0.1); link(3,4,0.5); link(4,5,0.5); link(2,3,0.1)
x0 = np.array([-1.6]+[-2.4]*5)
for K in [0.0, 0.5, 1.0, 2.0, 4.0]:
    o = simulate(C, x0, K=K, n_trials=4)
    print("K=",K, (o>0).any(axis=0).mean(0).round(2))
Cf, labels = load_template(); M, g = coarse(Cf, labels); M = np.minimum(M/M[0].max(),1.0)
x0 = -2.4 + 0.3*M[0]; x0[0] = -1.6
t = time.time()
for K in [0.3,0.5,0.7,1.0]:
    o = simulate(M, x0, K=K, n_trials=4, seed=5)
    print("10-group K", K, "->", int(((o>0).any(axis=0).mean(0)>=0.5).sum()), "of 10")
print("4 sims took", round(time.time()-t,1), "s")
