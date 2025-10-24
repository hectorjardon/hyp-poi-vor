# Plot Poisson-Voronoi on Euclidean space.

import numpy as np
import numpy.linalg as la   
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy.spatial import ConvexHull
from matplotlib.patches import Circle

# Generate Poisson point process in the Poincaré disk
def poisson_points_euc(intensity,radius): # intensity, radius of the disk of sampling in hyperbolic metric
    num_points = np.random.poisson(intensity * np.pi * radius**2) # number of points in the disk, Poison (intensity*area)
    radii_coordinates = np.random.uniform(low=0, high=1, size=num_points) # sample num_points from uniform in [0,1]
    radii = np.sqrt(radii_coordinates) * radius# sample radii by inverse transform sampling
    angles = np.random.uniform(low=0, high=2*np.pi, size=num_points) #angles sampled form uniform in [0,2pi]
    return np.array([radii*np.cos(angles),radii*np.sin(angles)]).T 
#output is points in Poincare disk with shape (num_points,2), for compatibility with Delaunay

def circumcenters_2d(triangles):

    A = triangles[:, 0, :] #Each A, B, C correspond to the vertices of the triangles
    B = triangles[:, 1, :]
    C = triangles[:, 2, :]
    
    # Compute determinants
    D = 2 * (A[:,0]*(B[:,1]-C[:,1]) + B[:,0]*(C[:,1]-A[:,1]) + C[:,0]*(A[:,1]-B[:,1]))
    
    # Circumcenter coordinates
    U = ((A[:,0]**2 + A[:,1]**2)*(B[:,1]-C[:,1]) +
         (B[:,0]**2 + B[:,1]**2)*(C[:,1]-A[:,1]) +
         (C[:,0]**2 + C[:,1]**2)*(A[:,1]-B[:,1])) / D

    V = ((A[:,0]**2 + A[:,1]**2)*(C[:,0]-B[:,0]) +
         (B[:,0]**2 + B[:,1]**2)*(A[:,0]-C[:,0]) +
         (C[:,0]**2 + C[:,1]**2)*(B[:,0]-A[:,0])) / D
    
    circumcenters = np.stack([U, V], axis=1)
    
    # Circumradius: distance from circumcenter to any vertex (here A)
    circumradii = np.linalg.norm(circumcenters - A, axis=1)
    
    return circumcenters, circumradii


def plot_euc_PVT(intensity,radius,tol):
    points = poisson_points_euc(intensity,radius + np.sqrt(-np.log(1-tol) / (intensity * np.pi))) # add margin to radius to avoid edge effects, ensures that with probability tol a point will be sampled in this windo
    tri = Delaunay(points) # Delaunay triangulation of the points
    triangles_idx = tri.simplices # indices of the triangles
    triangles_physical = points[triangles_idx] # shape (num_triangles,vertices,dimensions) where m is the number of triangles
    circumcenters, circumradii = circumcenters_2d(triangles_physical) # shape (num_triangles,2), (num_triangles,)

    mask = (triangles_idx[:,:,np.newaxis] == range(len(points))).any(axis=1) # mask to identify triangles containing each point
    hulls = [ConvexHull(circumcenters[mask[:,i],:]) for i in range(len(points)) if np.sum(mask[:,i])>=3] # only keep cells with at least 3 vertices
    patches = [Polygon(hull.points[hull.vertices],closed = True) for hull in hulls] # list of polygons for the cells
    colors = np.random.rand(len(points),4) # random colour for each cell
    collection = PatchCollection(patches, facecolor=colors, edgecolor='black', alpha=0.7)

    # Plot
    fig, ax = plt.subplots()
    ax.add_collection(collection)
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_aspect('equal','box')
    ax.axis('off')
    ax.set_title(f"Intensity = {intensity}")
    plt.show()
    return 

plot_euc_PVT(0.0001,5,0.999)