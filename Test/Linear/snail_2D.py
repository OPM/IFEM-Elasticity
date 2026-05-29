from sys import argv
from splipy.io import G2
import splipy.curve_factory as cf
import splipy.surface_factory as sf
import numpy as np

p = int(argv[1]) if len(argv) > 1 else 1

def helix1(t):
    return np.vstack([(t+1)*np.cos(t*2*np.pi), (t+1)*np.sin(t*2*np.pi)]).T
def helix2(t):
    return np.vstack([(t+2)*np.cos(t*2*np.pi), (t+2)*np.sin(t*2*np.pi)]).T

c1 = cf.fit(helix1, 0, 2)
c2 = cf.fit(helix2, 0, 2)
s1 = sf.edge_curves(c2, c1)

with G2("snail-p" + str(p) + ".g2") as output:
    output.write(s1.rebuild(p=p+1, n=(21*p,p+1)))
