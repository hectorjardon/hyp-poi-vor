import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
# We want to compute the hyperbolic distance between two points in the Poincaré disk model.
def hyperbolic_distance(p1,p2):
    p1,p2 = np.array(p1), np.array(p2)
    assert la.norm(p1) < 1 and la.norm(p2) < 1, "Points must be inside the unit disk"
    euc_dist = la.norm(p1 - p2)**2
    denom = (1 - la.norm(p1) **2) * (1 - la.norm(p2)**2)
    return np.arccosh(1 + 2 * euc_dist / denom)



def poisson_points(intensity,radius): #intensity, radius of the disk in hyperbolic metric
    num_points=np.random.poisson(intensity*2*np.pi*(np.cosh(radius)-1))
    radii_coordinates = np.random.uniform(low=0,high=1,size=num_points)
    radii = np.tanh(np.arccosh(1 + radii_coordinates*(np.cosh(radius)-1)))
    angles = np.random.uniform(low=0,high=2*np.pi,size=num_points)
    return np.array([[r*np.cos(a),r*np.sin(a)] for r,a in zip(radii,angles)])
#output is points in Poincare disk

    

def mesh(sqrt_points):
    radii = np.linspace(-0.5,0.5,sqrt_points)
    angles = np.linspace(0,2*np.pi,sqrt_points)
    R,A = np.meshgrid(radii,angles)
    X,Y = R*np.cos(A), R*np.sin(A)
    return X,Y

def colours(X,Y,points):
    colours = np.zeros_like(X)
    colour_vect = np.random.rand(points.shape[0])
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            dist = np.zeros_like(points[:,0],dtype=float)   
            for k in range(points.shape[0]):
                dist[k] = hyperbolic_distance([X[i,j],Y[i,j]],points[k])
            colours[i,j] = colour_vect[np.argmin(dist)]
    return colours
        
radius = 5
intensity = 0.1
w = np.sqrt((np.cosh(radius/10)-1)/(np.cosh(radius/10)+1))
x = np.linspace(-w,w,10)
y = np.linspace(-w,w,10)
X,Y = np.meshgrid(x,y)
points = poisson_points(intensity,radius)
Z = colours(X,Y,points)

plt.figure(figsize=(6,6))
plt.pcolormesh(X, Y, Z, cmap='viridis')  # 'viridis' is a nice colormap
plt.title('Poisson-Voronoi tessellation')
plt.gca().set_aspect('equal')  # keep axes equal
plt.show()


