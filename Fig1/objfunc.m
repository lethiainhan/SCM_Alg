function f = objfunc(x,y,FuncName)
% define a liste of 8 bivariate objectif functions for maximization testing

% x and y are coordinates and FuncName is the function name in string form

switch FuncName
    % function given from Bouman's book page 92 
    % => to show the diffences between GA, ICA and SCA (with zoom in)
    case 'Bouma'
        f = exp(-4*(x+y).^2)/((x-y).^2+1);
    % function given from https://en.wikipedia.org/wiki/Coordinate_descent
    % => to show the failure of GA, ICA and SCA
    case 'Wikip'
        f = -abs(x+y)-3*abs(y-x);
    % Rosenbrock function given from https://en.wikipedia.org/wiki/Rosenbrock_function
    % => to show the advantage of SCA over GA and ICA
    case 'Rosen'
        f = -(1-x)^2-100*(y-x^2)^2;
    % Sphere function given from https://www.sfu.ca/~ssurjano/spheref.html
    % => to show in the easy case where variables are independent
    case 'Spher'
        f = -x^2-y^2;
    % Sum of different powers function given from https://www.sfu.ca/~ssurjano/sumpow.html
    % => to show independent variables, but plat optimal zone (GA slow, others quick) 
    case 'Sumpo'
        f = -x^2-(abs(y))^3;
    % Booth function given from https://www.sfu.ca/~ssurjano/booth.html
    % => to show dependent variables, SCA is better than GA
    case 'Booth'
        f = -(x+2*y-7)^2-(2*x+y-5)^2; 
    % Matya function given from https://www.sfu.ca/~ssurjano/matya.html
    % => to show dependent variables, GA is better than SCA
    case 'Matya'
        f = -0.26*(x^2+y^2)+0.48*x*y;
    % Zakharov function given from https://www.sfu.ca/~ssurjano/zakharov.html
    case 'Zakha'
        t = 0.5*1*x + 0.5*2*y;
        f = -(x^2+y^2+t^2+t^4);
    % Mccormick function given from https://www.sfu.ca/~ssurjano/mccorm.html
    otherwise
        f = -sin(x+y) - (x-y)^2 + 1.5*x - 2.5*y - 1;
end

end