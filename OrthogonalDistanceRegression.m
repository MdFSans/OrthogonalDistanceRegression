function [vd,P] = OrthogonalDistanceRegression(vect)
    % Fitting a vector data 'vect' using Orthogonal Distance Regression. 
    % Minimize the ortogonal distance by adjusting both the values of
    % the dependent and independent variables.
    % Objective function to minimize: F = Sum(xi*sin(alpha) + yi*cos(alpha) + c )^2
    % There are two angle which achieve dF/dc = 0 and dF/dalpha = 0. One of
    % them maximizes the function and other minimizes it.
    
    % Initialize summations
    Sy = 0;
    Sx = 0;
    Sxy = 0;
    Sxx = 0;
    Syy = 0;
    ndata = size(vect,1);
    if ndata <= 1
        errordlg('It is not possible to adjust a single value')
    end
    for idata = 1 : ndata
        Sx = Sx + vect(idata,1);
        Sy = Sy + vect(idata,2);
        Sxy = Sxy + vect(idata,1)*vect(idata,2);
        Sxx = Sxx + vect(idata,1)*vect(idata,1);
        Syy = Syy + vect(idata,2)*vect(idata,2);
    end
    
    % Find the initial angle
    alpha1 = 1/2*atand(-2*(Sxy - 1/ndata*Sx*Sy)/(Sxx - Syy - 1/ndata*Sx*Sx + 1/ndata*Sy*Sy));
    if alpha1 > -1e-12 && alpha1 < 1e-12
        alpha1 = 0;
    end
    
    % Find the parameters of the function for both angles.
    alpha = zeros(2,1);
    c = zeros(2,1);
    F = zeros(2,1);
    for i = 1 : 2
        alphaAux = alpha1 + (i-1)*90;
        cAux = -1/ndata*Sx*sind(alphaAux) - 1/ndata*Sy*cosd(alphaAux);
        FAux = Sxx*(sind(alphaAux))^2 + Syy*(cosd(alphaAux))^2 + ndata*cAux^2 + ...
            2*(Sxy*sind(alphaAux)*cosd(alphaAux) + cAux*Sy*cosd(alphaAux) + cAux*Sx*sind(alphaAux));
        alpha(i,1) = alphaAux;
        c(i,1) = cAux;
        F(i,1) = FAux;
    end
    
    % Choose the angle which minimizes the function
    if F(1,1) < F(2,1) || F(1,1) == F(2,1)
        num = 1;
    elseif F(1,1) > F(2,1)
        num = 2;
    end
    % Find the slope and y-intercept.
    % Note that if alpha is 90º slope is infinite and n is x-intercept
    if alpha(num,1) ~= 90
        m = -sind(alpha(num,1))/cosd(alpha(num,1));
        n = -c(num,1)/cosd(alpha(num,1));
        vd = [1 m];
        P = [0 n];
    else
        m = 'infinite';
        n = -c(num,1)/sind(alpha(num,1));
        vd = [0 1];
        P = [n 0];
    end
end