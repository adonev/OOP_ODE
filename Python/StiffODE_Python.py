B = numpy.mat([0.2e0*x[0]*x[1]+0.748585e0*math.cos(t)+0.911111e0*math.sin(t)-0.2e0*math.sin(t)*math.cos(t),
               0.3e0*x[0]**2-0.1e0*x[1]**2-0.125142e1*math.sin(t)+0.288888e0*math.cos(t)-0.3e0*math.sin(t)**2+0.1e0*math.cos(t)**2])

A = numpy.mat([[0.222222e-1+0.888889e0*eig,0.628539e-1-0.314269e0*eig],[0.628536e-1-0.314269e0*eig,0.177777e0+0.111111e0*eig]])
