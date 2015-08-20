""" Test various methods for generating points on spheres """
import math
import random
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

def apbs_sphere(radius, count):
    """ This is the way APBS attempts to distribute points on a sphere more-or-less
    uniformly """
    frac = count/4.0
    ntheta = math.floor(math.sqrt(math.pi*frac)+0.5)
    dtheta = (math.pi)/(float(ntheta))
    nphimax = 2 * ntheta
    x = []
    y = []
    z = []
    nactual = 0
    for itheta in range(int(ntheta)):
        theta = dtheta * float(itheta)
        sintheta = math.sin(theta)
        costheta = math.cos(theta)
        nphi = math.floor((sintheta*nphimax)+0.5)
        nactual += nphi
        if nphi > 0:
            dphi = (2*math.pi)/float(nphi)
            for iphi in range (int(nphi)):
                phi = dphi*(float(iphi))
                sinphi = math.sin(phi)
                cosphi = math.cos(phi)
                xpts = radius * cosphi * sintheta
                ypts = radius * sinphi * sintheta
                zpts = radius * costheta
                x.append(xpts)
                y.append(ypts)
                z.append(zpts)
                nactual += 1
    return (x, y, z)

def random_sphere(radius, count):
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

def regular_sphere(radius, count):
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
        sintheta = math.sin(theta)
        costheta = math.cos(theta)
        M_phi = int(round(2*math.pi*sintheta/d_phi))
        for n in range(M_phi):
            phi = 2*math.pi*n/M_theta
            cosphi = math.cos(phi)
            sinphi = math.sin(phi)
            x_val = radius*(sintheta * cosphi)
            y_val = radius*(sintheta * sinphi)
            z_val = radius*(costheta)
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
    return (mean_d, var_d)

def scatter_plot3d(x, y, z, title=None):
    """ Scatter plot a set of points in 3D """
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
    count = 5000
    test_radius = 0.001 * radius

    print "Testing the random code..."
    (sx, sy, sz) = random_sphere(radius, count)
    scatter_plot3d(sx, sy, sz, "Random")
    print "Density mean = %g, variance = %g" % density_variance(sx, sy, sz, test_radius)

    print "Testing the new uniform code"
    (rx, ry, rz) = regular_sphere(radius, count)
    scatter_plot3d(rx, ry, rz, "New")
    print "Density mean = %g, variance = %g" % density_variance(rx, ry, rz, test_radius)

    print "Testing the APBS uniform code"
    (ax, ay, az) = apbs_sphere(radius, count)
    scatter_plot3d(ax, ay, az, "APBS")
    print "Density mean = %g, variance = %g" % density_variance(ax, ay, az, test_radius)
