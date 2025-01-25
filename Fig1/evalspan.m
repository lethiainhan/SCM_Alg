function [xmin,xmax,ymin,ymax] = evalspan(FuncName)
% define a liste of 8 bivariate objectif functions for maximization testing

% x and y are coordinates and FuncName is the function name in string form

switch FuncName
    % function given from Bouman's book page 92
    case 'Bouma'
        xmin = -2; xmax = 2; ymin = -2; ymax = 2;
    % function given from https://en.wikipedia.org/wiki/Coordinate_descent
    case 'Wikip'
        xmin = -3; xmax = 3; ymin = -3; ymax = 3;
    % Rosenbrock function given from https://en.wikipedia.org/wiki/Rosenbrock_function
    case 'Rosen'
        xmin = -5; xmax = 5; ymin = -5; ymax = 5;
    % Sphere function given from https://www.sfu.ca/~ssurjano/spheref.html
    case 'Spher'
        xmin = -5; xmax = 5; ymin = -5; ymax = 5;
    % Sum of different powers function given from https://www.sfu.ca/~ssurjano/sumpow.html
    case 'Sumpo'
        xmin = -1; xmax = 1; ymin = -1; ymax = 1;
    % Booth function given from https://www.sfu.ca/~ssurjano/booth.html
    case 'Booth'
        xmin = -10; xmax = 10; ymin = -10; ymax = 10;
    % Booth function given from https://www.sfu.ca/~ssurjano/matya.html
    case 'Matya'
        xmin = -10; xmax = 10; ymin = -10; ymax = 10;
    % Zakharov function given from https://www.sfu.ca/~ssurjano/zakharov.html
    case 'Zakha'
        xmin = -10; xmax = 10; ymin = -10; ymax = 10;
    % Mccormick function given from https://www.sfu.ca/~ssurjano/mccorm.html
    otherwise
        xmin = -4; xmax = 4; ymin = -4; ymax = 4;
end

end