def scatter_regular():
    import math
    N_count = 0
    count = 5000
    r = 1.0
    a = 4*math.pi*math.pow(r,2)/count
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
            x_val = r*(math.sin(theta) * math.cos(phi))
            y_val = r*(math.sin(theta) * math.sin(phi))
            z_val = r*(math.cos(theta))
            x.append(x_val)
            y.append(y_val)
            z.append(z_val)
            N_count += 1
    print len(x)
    scatter(x, y, z)
    return(x, y, z)
if __name__ == 'main':
    scatter_regular()