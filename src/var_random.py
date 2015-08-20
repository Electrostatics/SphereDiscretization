def var_random():
    (x, y, z) = scatter_random()
    rA = .001 * r
    d = [0.0] * count
    rA2 = math.pow(rA,2)
    for i in range(0, count):
        for j in range(0, count):
            if rA2 >= (x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + (z[i] - z[j]) ** 2:
                d[i] += 1.0

    print d
    #print [x for x in d if x > 1.0]
    mean_d = sum(d)
    mean_d /= float(len(d))

    var_d = sum((x-mean_d)**2 for x in d)

    var_d /= float(len(d))
    print(mean_d, var_d)
if __name__ == 'main':
    var_random()