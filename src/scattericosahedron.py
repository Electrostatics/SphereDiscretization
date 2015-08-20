""" Test various methods for generating points on spheres """
import math
import numpy as np
import random
import argparse
import os
import sys
import math
import time
import sys
from sys import argv
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D


def scatter_icosahedron(filename):
    """ This function takes the points from the txt file created by subdiver.py to use for graphing """
    #script, filename = argv

    x, y, z = np.genfromtxt(filename, unpack=True)
    #print x
    return(x, y, z)

def density_variance(x, y, z, test_radius):
    """ This attempts to estimate the heterogeneity of points on the surface of
    a sphere """
    count = len(x)
    d = [0.0] * count
    tr2 = test_radius * test_radius
    # Count the number of points within/outside a certain distance of the
    # current point
    for i in range(0, count):
        for j in range(i+1, count):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            dz = z[i] - z[j]
            if dx*dx + dy*dy + dz*dz <= tr2:
                d[i] += 1.0
    mean_d = sum(d)
    mean_d = mean_d/float(len(d))
    var_d = sum((x-mean_d)*(x-mean_d) for x in d)
    var_d = var_d/float(len(d))
    return (mean_d, var_d)

def scatter_plot3d(x, y, z, title=None):
    """Scatter plot a set of points in 3D"""
    pyplot.ioff()
    fig = pyplot.figure()
    axes = fig.gca(projection='3d')
    axes.scatter(xs=x, ys=y, zs=z)
    if title:
        pyplot.title(title)
        pyplot.show()


if __name__ == '__main__':
    matplotlib.use('macosx')
    pyplot.close("all")

    # These are the user-customizable settings
    radius = 1.0
    test_radius = .05 * radius
    filename = sys.argv[1]
    print "Testing the icosahedron code"
    (mx, my, mz) = scatter_icosahedron(filename)
    scatter_plot3d(mx, my, mz, "Icosahedron")
    print "Density mean = %g, variance = %g" %density_variance(mx, my, mz, test_radius)

