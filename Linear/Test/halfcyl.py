import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory
from splipy.io import G2
from math import pi

R    = 0.1
L    = 2.0
arch = CurveFactory.circle_segment(pi,R,[3,0,0],normal=(1,0,0))
surf = SurfaceFactory.extrude(arch,amount=(L,0,0))

with G2('halfcyl.g2') as output:
    output.write(surf)
