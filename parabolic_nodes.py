
#!/usr/bin/python

import sys
import math


def f(Rmin, Rmax, R0, linear_slope):
    i0 = math.floor(2*(R0-Rmin)/linear_slope)
    nodes = (2*(Rmax - Rmin)/linear_slope + i0)/2
    print(f"To have a slope of <= {linear_slope} at R={R0} you need >={nodes} nodes")

f(float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]))