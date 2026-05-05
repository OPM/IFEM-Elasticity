import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory
from splipy.io import G2
from math import sqrt, sin, pi

beta = 0.1
R = 2540
l = 508
w = 2*R*sin(beta)
y0 = -sqrt(R*R-0.25*w*w)
arch = CurveFactory.circle_segment(2*beta,R,[0,0,1])
arch.rotate(pi/2-beta)
arch.translate([0,y0,0])
shell = SurfaceFactory.extrude(arch,[0,0,l])
print(shell)

with G2('shallow_arch.g2') as output:
    output.write(shell)
