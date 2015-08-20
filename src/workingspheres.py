import math
import random
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

def scatter_random(radius, count):
    """ This function generates <count> points on a sphere of radius
    <radius> """
    x = []
    y = []
    z = []
    for a in range(count):
        z_term = random.uniform(-radius, radius)
        phi = random.uniform(0, 2 * math.pi)
        x_term = math.sqrt((radius * radius - z_term * z_term)) * math.cos(phi)
        y_term = math.sqrt((radius * radius - z_term * z_term)) * math.sin(phi)
        x.append(x_term)
        y.append(y_term)
        z.append(z_term)
    return  (x, y, z)

def scatter_regular(radius, count):
    """ This function generates approximately <count> points on a sphere of
    radius <radius """
    N_count = 0
    a = 4*math.pi*radius*radius/float(count)
    d = math.sqrt(a)
    M_theta = int(round(math.pi/d))
    d_theta = math.pi/M_theta
    d_phi = a/d_theta
    x = []
    y = []
    z = []
    for m in range(M_theta):
        theta = math.pi*(m + .5)/M_theta
        M_phi = int(round(2*math.pi*math.sin(theta)/d_phi))
        for n in range(M_phi):
            phi = 2*math.pi*n/M_theta
            x_val = radius*(math.sin(theta) * math.cos(phi))
            y_val = radius*(math.sin(theta) * math.sin(phi))
            z_val = radius*(math.cos(theta))
            x.append(x_val)
            y.append(y_val)
            z.append(z_val)
            N_count += 1
    return (x, y, z)

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
    print(mean_d, var_d)

if __name__ == '__main__':
    matplotlib.use('macosx')
    pyplot.close("all")

    # These are the user-customizable settings
    radius = 1.0
    count = 5000
    test_radius = 0.001 * radius

    print "Testing the random code"
    (sx, sy, sz) = scatter_random(radius, count)
    pyplot.ioff()
    pyplot.figure()
    pyplot.title("Random")
    pyplot.scatter(sx, sy)
    pyplot.show()
    density_variance(sx, sy, sz, test_radius)

    print "Testing the regular code"
    (rx, ry, rz) = scatter_regular(radius, count)
    pyplot.ioff()
    pyplot.figure()
    pyplot.title("Regular")
    pyplot.scatter(rx, ry)
    pyplot.show()
    density_variance(rx, ry, rz, test_radius)
    print rz[80:100]
