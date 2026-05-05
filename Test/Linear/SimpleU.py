import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory
import splipy.volume_factory as VolumeFactory
from splipy.io import G2

L = 5.0
W = 0.6
H = 0.4
t = 0.1
patch = [None] * 3
lines = [None] * 6
lines[0] = CurveFactory.line([0,0,H],[0,0,0])
lines[1] = CurveFactory.line([0,0,0],[0,W,0])
lines[2] = CurveFactory.line([0,W,0],[0,W,H])
lines[3] = CurveFactory.line([0,t,H],[0,t,t])
lines[4] = CurveFactory.line([0,t,t],[0,W-t,t])
lines[5] = CurveFactory.line([0,W-t,t],[0,W-t,H])
for i in range(3):
    surface = SurfaceFactory.loft(lines[i],lines[i+3])
    patch[i] = VolumeFactory.extrude(surface,[L,0,0])

with G2('Uprofil.g2') as output:
    output.write(patch)
