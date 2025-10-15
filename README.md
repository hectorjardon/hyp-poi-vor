This projects provides Python functions to experiment with Poisson--Voronoi tesselations and their low intensity limits. 
Properties of the low intensity tesselations are closely related to the geometry of the space as has been notably manifested in the work of Fraczyk-Mellick-Wilkens, or in my work with Mellick-Poulin-Wrobel.

The file plot_PVT.py contains a function plot_PVT which plots the Poisson--Voronoi tesselation in the Poincare disk model of hyperbolic geometry. It produces the plot using a mesh, so it is precise, but costly.

The file klein_PVT.py contains plot_Klein_PVT, which does the same job on the Beltrami--Klein model of hyperbolic geometry. 
The plot is much faster as it relies on the Delaunay function of scipy, but works only for visualization purposes. 
In this file also a function plot_euc_PVT is present, plotting the Poisson--Voronoi tesselation in Euclidean plane for comparison.
