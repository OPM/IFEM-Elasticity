import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory
import splipy.volume_factory as VolumeFactory
from splipy.io import G2
from math import pi
from sys import argv

r1 = 5
r2 = 7
H  = 1

inner = CurveFactory.circle_segment(pi,r1,[0,0,1])
outer = CurveFactory.circle_segment(pi,r2,[0,0,1])
surf  = SurfaceFactory.edge_curves(outer,inner)
vol   = None

if len(argv) > 1:
    if argv[1] == "3D":
        vol = VolumeFactory.extrude(surf,amount=(0,0,H))
    g2file = "annulus" + argv[1] + ".g2"
else:
    g2file = "annulus.g2"

with G2(g2file) as output:
    if vol:
        output.write(vol)
    else:
        output.write(surf)
