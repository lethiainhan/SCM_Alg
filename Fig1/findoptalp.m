function alpopt = findoptalp(x0,y0,dfx0,dfy0,FuncName)
% find optimal value of alpha in line search using root of derivative function wrt alpha

% define constants
diffdf = dfx0 - dfy0; sumdf = dfx0 + dfy0; diffxy0 = x0 - y0; sumxy0 = x0 + y0;

% define the derivative function wrt alpha
switch FuncName
    % function given from Bouman's book page 92
    case 'Bouma'
        devfalp = @(alp) diffdf*(diffdf*alp + diffxy0) + 4*sumdf*((diffdf*alp + diffxy0)^2 + 1)*(sumdf*alp + sumxy0);
    % function given from https://en.wikipedia.org/wiki/Coordinate_descent
    case 'Wikip'
        devfalp = @(alp) -sumdf*(sumdf*alp + sumxy0)/abs(sumdf*alp + sumxy0 +eps) ...
                   - 3*diffdf*(diffdf*alp + diffxy0)/abs(diffdf*alp + diffxy0 +eps);
    % Rosenbrock function given from https://en.wikipedia.org/wiki/Rosenbrock_function
    case 'Rosen'
        devfalp = @(alp) 2*dfx0*(1 - x0 - dfx0*alp) - 200*(dfy0 - 2*dfx0*(x0 + dfx0*alp))...
                   * (y0 + dfy0*alp - (x0 + dfx0*alp)^2);
    % Sphere function given from https://www.sfu.ca/~ssurjano/spheref.html
    case 'Spher'
        devfalp = @(alp) -2*dfx0*(x0 + dfx0*alp) - 2*dfy0*(y0 + dfy0*alp);
    % Sum of different powers function given from https://www.sfu.ca/~ssurjano/sumpow.html
    case 'Sumpo'
        devfalp = @(alp) -2*dfx0*(x0 + dfx0*alp) - 3*dfy0*(y0 + dfy0*alp)*abs(y0 + dfy0*alp);
    % Booth function given from https://www.sfu.ca/~ssurjano/booth.html
    case 'Booth'
        devfalp = @(alp) -10*dfx0^2*alp - 2*dfy0*(-19+4*x0+5*y0+5*dfy0*alp)...
                    - 2*dfx0*(-17+5*x0+4*y0+8*dfy0*alp); 
    % Booth function given from https://www.sfu.ca/~ssurjano/matya.html
    case 'Matya'
        devfalp = @(alp) -0.52*dfx0*(x0 + dfx0*alp) + 0.48*dfy0*(x0 + dfx0*alp)...
                    + 0.48*dfx0*(y0 + dfy0*alp) - 0.52*dfy0*(y0 + dfy0*alp);
    % Zakharov function given from https://www.sfu.ca/~ssurjano/zakharov.html
    case 'Zakha'
        devfalp = @(alp) -2*dfx0*(x0 + dfx0*alp) - 2*dfy0*(y0 + dfy0*alp)...
                    - ((dfx0 + 2*dfy0)*(x0 + 2*y0 + (dfx0 + 2*dfy0)*alp))/2 ...
                    - ((dfx0 + 2*dfy0)*(x0 + 2*y0 + (dfx0 + 2*dfy0)*alp)^3)/4;
    % Mccormick function given from https://www.sfu.ca/~ssurjano/mccorm.html
    otherwise
        devfalp = @(alp) 1.5*dfx0 - 2.5*dfy0 - 2*(dfx0 - dfy0)*(diffxy0 + diffdf*alp)...
                    - sumdf*cos(sumxy0 + sumdf*alp);
end

% Find the root of the simplied number function
alpopt = fzero(devfalp,1);

end