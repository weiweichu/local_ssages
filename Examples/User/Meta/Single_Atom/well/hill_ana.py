import numpy as np 
import matplotlib.pyplot as plt
def min_dist(x):
     return np.mod(x+1.75,3.5)-1.75
def sum_hills(hills, first, last, points, lbnds, ubnds): 
     # Get number of CVs.
    ncv = int((hills.shape[1] - 2)/2)
    centers = hills[first:last,1:ncv+1]
    sigmas  = hills[first:last,ncv+1:2*ncv+1]
    height  = hills[first:last,-1]

        # Grid data and sum hills.
    X = np.meshgrid(*[np.linspace(i,j,n) for i,j,n in zip(lbnds,ubnds,points)])
    X = np.array([x.flatten() for x in X]).T
    # Evaluate hills.
    Z = np.zeros((X.shape[0],))
    for i in range(len(Z)):
        Z[i] = np.sum(height*np.prod([np.exp(-min_dist(x-c)**2/(2.*s**2)) \
                for x,c,s in zip(X[i,:],centers.T,sigmas.T)],axis=0))

    return X, Z
dt = 0.001 # Timestep in fs. 
hfreq = 500 # Frequency of hill drops (in iterations)
times = np.array([1000, 2000, 4000]) # Times at which to plot FES (ps).

v = np.arange(-50, 5, 2) # Contours to plot.
for t in times:
    fig = plt.figure(figsize=(10, 8.2))
    frame = int(t/(hfreq*dt))
    hills = np.loadtxt("hills.out".format(frame), skiprows=1)

        # Reshape and sum hills.
    x, z = sum_hills(hills, 0, frame, [50, 50], [-1.5,-1.5], [2.0, 2.0])
    xg = x[:,0].reshape((50, 50))
    yg = x[:,1].reshape((50, 50))
    zg = -z.reshape((50, 50))
    zg = zg - np.max(zg)
        
    # Plot data.
    plt.title("Metadynamics at t = {} ns".format(t/1000.))
    plt.contour(xg, yg, zg, v, linewidths=0.5, colors="k")
    #plt.contourf(xg, yg, zg, cmap=plt.cm.plasma)
    plt.contourf(xg, yg, zg, v, cmap=plt.cm.plasma)
    cb = plt.colorbar()
    cb.set_label("G(kJ/mol)")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.xlim((-1.5, 2.0))
    plt.ylim((-1.5,2.0))
    plt.axis("equal")
    plt.savefig("meta_{:.1f}ns.png".format(t/1000.))
