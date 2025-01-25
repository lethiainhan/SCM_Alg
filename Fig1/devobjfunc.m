function [dfx0,dfy0] = devobjfunc(x0,y0,FuncName)
% compute partial derivatives of the objectif function objfunc(x,y) at (x0,y0)

% ================ if numerical computation used ================
% %
% % define steps
% delx = 1e-10; dely = 1e-10;
% 
% % partial derivative with respect to x
% dfx0 = (objfunc(x0+delx,y0,FuncName) - objfunc(x0,y0,FuncName))/delx;
% 
% % partial derivative with respect to y
% dfy0 = (objfunc(x0,y0+dely,FuncName) - objfunc(x0,y0,FuncName))/dely;

% ================ if exact computation used ================

switch FuncName
    % function given from Bouman's book page 92
    case 'Bouma'
        dfx0 = -2*(5*x0 + 4*x0^3 + 3*y0 - 4*x0^2*y0 - 4*x0*y0^2 + 4*y0^3)/(exp(4*(x0 + y0)^2)*(1 + x0^2 - 2*x0*y0 + y0^2)^2);
        dfy0 = -2*(3*x0 + 4*x0^3 + 5*y0 - 4*x0^2*y0 - 4*x0*y0^2 + 4*y0^3)/(exp(4*(x0 + y0)^2)*(1 + x0^2 - 2*x0*y0 + y0^2)^2);
    % function given from https://en.wikipedia.org/wiki/Coordinate_descent
    case 'Wikip'
        dfx0 = 3*(-x0 + y0)/abs(x0 - y0) - (x0 + y0)/abs(x0 + y0);
        dfy0 = 3*(x0 - y0)/abs(x0 - y0) - (x0 + y0)/abs(x0 + y0);
    % Rosenbrock function given from https://en.wikipedia.org/wiki/Rosenbrock_function
    case 'Rosen'
        dfx0 = 2 - 400*x0^3 + x0*(-2 + 400*y0);
        dfy0 = -200*(-x0^2 + y0);
    % Sphere function given from https://www.sfu.ca/~ssurjano/spheref.html
    case 'Spher'
        dfx0 = -2*x0;
        dfy0 = -2*y0;
    % Sum of different powers function given from https://www.sfu.ca/~ssurjano/sumpow.html
    case 'Sumpo'
        dfx0 = -2*x0;
        dfy0 = -3*y0*abs(y0);
    % Booth function given from https://www.sfu.ca/~ssurjano/booth.html
    case 'Booth'
        dfx0 = -10*x0 - 8*y0 + 34;
        dfy0 = -8*x0 - 10*y0 + 38;
    % Booth function given from https://www.sfu.ca/~ssurjano/matya.html
    case 'Matya'
        dfx0 = -0.52*x0 + 0.48*y0;
        dfy0 = 0.48*x0 - 0.52*y0;
    % Zakharov function given from https://www.sfu.ca/~ssurjano/zakharov.html
    case 'Zakha'
        dfx0 = -5*x0/2 - y0 - (x0 + 2*y0)^3/4;
        dfy0 = 2*(-x0/2 - 2*y0 - (x0 + 2*y0)^3/4);
    % Mccormick function given from https://www.sfu.ca/~ssurjano/mccorm.html
    otherwise
        dfx0 = -2*(-0.75 + x0 - y0 + 0.5*cos(x0 + y0));
        dfy0 = -2.5 + 2*x0 - 2*y0 - cos(x0 + y0);
end
end