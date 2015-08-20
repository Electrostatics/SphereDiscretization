def scatter_random():
    a = 0
    count = 5000
    r = 1.0
    x = []
    y = []
    z = []
    import math
    import random
    for a in range(count):
        z_term = random.uniform(-r, r)
        phi = random.uniform(0, 2 * math.pi)
        x_term = math.sqrt((r * r - z_term * z_term)) * math.cos(phi)
        y_term = math.sqrt((r * r - z_term * z_term)) * math.sin(phi)
        x.append(x_term)
        y.append(y_term)
        z.append(z_term)
    #print len(x)
    scatter(x, y, z)
    return(x, y, z)
if __name__ == 'main':
    scatter_random()