function [F,X,Y] = objfuncmesh(n,FuncName)
% define the rectangular grid and values of 8 bivariate objectif functions
% n is the number of points in grid
% example
% n = 101; [F,X,Y] = objfuncmesh(n,'Bouma'); 
% figure; mesh(X,Y,F); figure; contour(X,Y,F,50);

% x and y are coordinates and FuncName is the function name in string form

switch FuncName
    % function given from Bouman's book page 92
    case 'Bouma'
        % get x-y span 
        [xmin,xmax,ymin,ymax] = evalspan(FuncName);
        % define the rectangular grid
        [X,Y] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));
        % compute mesh of objectif junction
        F = exp(-4*(X+Y).^2)./((X-Y).^2+1);
    % function given from https://en.wikipedia.org/wiki/Coordinate_descent
    case 'Wikip'
        % get x-y span 
        [xmin,xmax,ymin,ymax] = evalspan(FuncName);
        % define the rectangular grid
        [X,Y] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));
        % compute mesh of objectif junction
        F = -abs(X+Y)-3*abs(Y-X);
    % Rosenbrock function given from https://en.wikipedia.org/wiki/Rosenbrock_function
    case 'Rosen'
        % get x-y span 
        [xmin,xmax,ymin,ymax] = evalspan(FuncName);
        % define the rectangular grid
        [X,Y] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));
        % compute mesh of objectif junction
        F = -(1-X).^2-100*(Y-X.^2).^2;
    % Sphere function given from https://www.sfu.ca/~ssurjano/spheref.html
    case 'Spher'
        % get x-y span 
        [xmin,xmax,ymin,ymax] = evalspan(FuncName);
        % define the rectangular grid
        [X,Y] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));
        % compute mesh of objectif junction        
        F = -X.^2-Y.^2;
    % Sum of different powers function given from https://www.sfu.ca/~ssurjano/sumpow.html
    case 'Sumpo'
        % get x-y span 
        [xmin,xmax,ymin,ymax] = evalspan(FuncName);
        % define the rectangular grid
        [X,Y] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));
        % compute mesh of objectif junction        
        F = -X.^2-(abs(Y)).^3;
    % Booth function given from https://www.sfu.ca/~ssurjano/booth.html
    case 'Booth'
        % get x-y span 
        [xmin,xmax,ymin,ymax] = evalspan(FuncName);
        % define the rectangular grid
        [X,Y] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));
        % compute mesh of objectif junction        
        F = -(X+2*Y-7).^2-(2*X+Y-5).^2; 
    % Booth function given from https://www.sfu.ca/~ssurjano/matya.html
    case 'Matya'
        % get x-y span 
        [xmin,xmax,ymin,ymax] = evalspan(FuncName);
        % define the rectangular grid
        [X,Y] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));
        % compute mesh of objectif junction        
        F = -0.26*(X.^2+Y.^2)+0.48*X.*Y;
    % Zakharov function given from https://www.sfu.ca/~ssurjano/zakharov.html
    case 'Zakha'
        % get x-y span 
        [xmin,xmax,ymin,ymax] = evalspan(FuncName);
        % define the rectangular grid
        [X,Y] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));
        % compute mesh of objectif junction 
        T = 0.5*1*X + 0.5*2*Y;
        F = -(X.^2+Y.^2+T.^2+T.^4);
    % Mccormick function given from https://www.sfu.ca/~ssurjano/mccorm.html
    otherwise
        % get x-y span 
        [xmin,xmax,ymin,ymax] = evalspan(FuncName);
        % define the rectangular grid
        [X,Y] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));
        % compute mesh of objectif junction        
        F = -sin(X+Y) - (X-Y).^2 + 1.5*X - 2.5*Y - 1;
end

end