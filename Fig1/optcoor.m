function [xopt,yopt] = optcoor(FuncName)
% define a liste of optimal coordinates of 8 bivariate objectif functions

% xopt and yopt are optimal coordinates of functions named by FuncName

switch FuncName
    % function given from Bouman's book page 92
    case 'Bouma'
        xopt = 0; yopt = 0;
    % function given from https://en.wikipedia.org/wiki/Coordinate_descent
    case 'Wikip'
        xopt = 0; yopt = 0;
    % Rosenbrock function given from https://en.wikipedia.org/wiki/Rosenbrock_function
    case 'Rosen'
        xopt = 1; yopt = 1;
    % Sphere function given from https://www.sfu.ca/~ssurjano/spheref.html
    case 'Spher'
        xopt = 0; yopt = 0;
    % Sum of different powers function given from https://www.sfu.ca/~ssurjano/sumpow.html
    case 'Sumpo'
        xopt = 0; yopt = 0;
    % Booth function given from https://www.sfu.ca/~ssurjano/booth.html
    case 'Booth'
        xopt = 1; yopt = 3;
    % Booth function given from https://www.sfu.ca/~ssurjano/matya.html
    case 'Matya'
        xopt = 0; yopt = 0;
    % Zakharov function given from https://www.sfu.ca/~ssurjano/zakharov.html
    case 'Zakha'
        xopt = 0; yopt = 0;
    % Mccormick function given from https://www.sfu.ca/~ssurjano/mccorm.html
    otherwise
        xopt = -0.54719; yopt = -1.54719;
end

end