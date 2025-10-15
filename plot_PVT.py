# This code implements a function plot_PVT that generates and visualizes a Poisson-Voronoi tessellation in the Poincaré disk model of hyperbolic geometry.
# The plot is constructed by a grid method.


import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

# We want to compute the hyperbolic distance between two points of the unit disk in the Poincaré disk model.
# Designed for input of arrays (2,n,n)
def hyperbolic_distance(p1,p2):
    assert (la.norm(p1,axis=0) < 1).all() and (la.norm(p2,axis=0) < 1).all(), "Points must be inside the unit disk"
    euc_dist = la.norm(p1 - p2,axis=0)**2
    denom = (1 - la.norm(p1,axis=0) **2) * (1 - la.norm(p2,axis=0)**2)
    return np.arccosh(1 + 2 * euc_dist / denom)
# Maybe useful to use the arctanh with complex numbers?

# Generate Poisson point process in the Poincaré disk
def poisson_points(intensity,radius): # intensity, radius of the disk of sampling in hyperbolic metric
    num_points = np.random.poisson(intensity * 2*np.pi * (np.cosh(radius)-1)) # number of points in the disk, Poison (intensity*area)
    radii_coordinates = np.random.uniform(low=0, high=1, size=num_points) # sample num_points from uniform in [0,1]
    radii = np.tanh(np.arccosh(1 + radii_coordinates*(np.cosh(radius)-1))/2) # sample radii by inverse transform sampling
    angles = np.random.uniform(low=0, high=2*np.pi, size=num_points) #angles sampled form uniform in [0,2pi]
    return np.array([radii*np.cos(angles),radii*np.sin(angles)]) 
#output is points in Poincare disk with shape (2,num_points)
    
# The input are two objects of shape (2,n,n) and (2,m). Want to compute distances for each (2,n,n) to every (2,m)
def dist_to_centers(grid,centers):
    centers_broad = np.broadcast_to(centers[:,:,np.newaxis,np.newaxis],(2,centers.shape[1],grid.shape[1],grid.shape[2]))
    grid_broad = np.broadcast_to(grid[:,np.newaxis,:,:],(2,centers.shape[1],grid.shape[1],grid.shape[2]))
    return hyperbolic_distance(grid_broad,centers_broad)
# Result is of shape (m,n,n) where m is the number of centers and n,n the grid shape and gives the distance to the corresponding center

# The input are two objects of shape (2,n,n) and (2,m). Output is (n,n) with the index of the closest center
def voronoi_center(grid,centers):
    return np.argmin(dist_to_centers(grid,centers),axis=0)

# Colour the Voronoi cells
def colours(grid,centers):
    colour_list = np.random.rand(centers.shape[1])
    id_center = voronoi_center(grid,centers)
    return colour_list[id_center]
# Output is (n,n) array with the colour of the closest center
        
# Plotting function

def plot_PVT(intensity,radius,density,tol):
    r = np.sqrt(np.linspace(0, np.tanh(radius/2)**2, density)) #radius in poincare disk
    a = np.linspace(0, 2*np.pi, density) #radius in poincare disk
    X = r[:,np.newaxis]*np.cos(a)[np.newaxis,:] # shape (n,n)
    Y = r[:,np.newaxis]*np.sin(a)[np.newaxis,:]
    mesh = np.stack((X,Y),axis=0) # shape (2,n,n)
    points = poisson_points(intensity,radius + np.arccosh(1 - np.log(1 - tol)/(2*np.pi*intensity))) # add margin to radius to avoid edge effects, ensures that with probability tol a point will be sampled in this window
    Z = voronoi_center(mesh,points)
    plt.figure(figsize=(6,6))
    plt.pcolormesh(X, Y, Z, cmap='gist_ncar')  # 'viridis' is a nice colormap
    plt.title('Poisson-Voronoi tessellation')
    plt.gca().set_aspect('equal')  # keep axes equal
    plt.show()
    return 

plot_PVT(1,5,300,0.9)



