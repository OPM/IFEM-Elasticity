import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory
import splipy.volume_factory as VolumeFactory
from splipy.io import G2
from math import pi

R    = 10
H    = 100
line = CurveFactory.line([0,0,0],[0,0,0])
arch = CurveFactory.circle_segment(0.5*pi,R,[0,0,1])
surf = SurfaceFactory.edge_curves(arch,line)
pch1 = VolumeFactory.extrude(surf,amount=(0,0,H))
pch2 = pch1.clone().rotate(0.5*pi)
pch3 = pch1.clone().rotate(1.0*pi)
pch4 = pch1.clone().rotate(1.5*pi)

with G2('bolt.g2') as output:
    output.write(pch1)
    output.write(pch2)
    output.write(pch3)
    output.write(pch4)
