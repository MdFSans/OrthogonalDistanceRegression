# OrthogonalDistanceRegression
Fitting a vector data 'vect' using Orthogonal Distance Regression (or Total least squares).

Minimize the ortogonal distance by adjusting both the values of the dependent and independent variables.

Objective function to minimize: F = Sum(xi·sin(alpha) + yi·cos(alpha) + c )^2
There are two angle which achieve dF/dc = 0 and dF/dalpha = 0. One ofthem maximizes the function and other minimizes it.
