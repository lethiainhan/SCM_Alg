close all

% ==== parameters definition ==== %
% define index for 2d and 3d plot: 1 = 2d and 0 = 3d ==== %
idx2d = 1;

% define the number of iterations used in algorithms ==== %
niter = 25;

% define an objectif function among 'Bouma', 'Wikip', 'Rosen', 'Spher',
% 'Sumpo', 'Booth', 'Matya', 'Zakha', 'Mccor'
FuncName = 'Bouma';

% get optimal coordinates
[xopt,yopt] = optcoor(FuncName);

% get x-y span
[xmin,xmax,ymin,ymax] = evalspan(FuncName);

% define randomly initial point using uniform distribution (see help rand in matlab)
% xinit = xmin + (xmax-xmin).*rand; yinit= ymin + (ymax-ymin).*rand;
xinit = -1.611; yinit = 1.294;


% define the rectangular grid and compute the associated values
n = 101; [F,X,Y] = objfuncmesh(n,FuncName);

% mesh plot
nf=900; figure(nf); mesh(X,Y,F); 
xlabel('x','FontSize',12), ylabel('y','FontSize',12)
set(gca,'FontSize',12)

% ==== optimization using gradient ascent with line search ==== %

% define the maximization algorithm's name
OptAlg = 'GAFR'; % (i.e., GAFR algorithm)

% preset (x0,y0)
x0 = xinit; y0 = yinit;

% plot the objective function, initial point and optimal point
nf=1000; FigContMesh(nf,X,Y,F,x0,y0,xopt,yopt,idx2d,FuncName);

% preallocation of residues and first residue
dGA = zeros(1,niter+1); 
dGA(1) = norm([xinit-xopt,yinit-yopt,objfunc(xinit,yinit,FuncName)-objfunc(xopt,yopt,FuncName)]);

% run gradient ascent with line search and plot
for i = 1:niter
    % compute partial derivative at (x0,y0) (eq. 5.23, page 97 of bouman's book)
    [dfx0,dfy0] = devobjfunc(x0,y0,FuncName);
    % find the optimal step size (eq. 5.24, page 97 of bouman's book)
    alpopt = findoptalp(x0,y0,dfx0,dfy0,FuncName);
    % update value (eq. 5.25, page 97 of bouman's book)
    x = x0 + alpopt*dfx0; 
    y = y0 + alpopt*dfy0;
    % compute the residue
    d = norm([x-xopt,y-yopt,objfunc(x,y,FuncName)-objfunc(xopt,yopt,FuncName)]); % dGA(i+1) = d;
    dGA(i+1) = objfunc(xopt,yopt,FuncName) - objfunc(x,y,FuncName);
    % plot the current updated point
    FigCurrPoinPlot(nf,x0,y0,x,y,idx2d,i,OptAlg,d,FuncName)
    % reset (x0,y0)
    x0 = x; y0 = y;
end

% ==== optimization using iterative coordinate ascent with line search ==== %

% define the maximization algorithm's name
OptAlg = 'ICMC'; % (i.e., ICMC algorithm)

% preset (x0,y0)
x0 = xinit; y0 = yinit;

% plot the objective function, initial point and optimal point
nf=2000; FigContMesh(nf,X,Y,F,x0,y0,xopt,yopt,idx2d,FuncName);

% preallocation of residues and first residue
dICA = zeros(1,niter+1); 
dICA(1) = norm([xinit-xopt,yinit-yopt,objfunc(xinit,yinit,FuncName)-objfunc(xopt,yopt,FuncName)]);

% define index for switch between x and y coordianate: 1 = x and 0 = y
idx = 1;

% run iterative coordinate ascent with line search and plot
for i = 1:niter
    if idx == 1
        % find the optimal step size following x coordinate
        alpoptx = findoptalp(x0,y0,1,0,FuncName);
        % update value
        x = x0 + alpoptx; 
        y = y0;
        % compute the residue
        d = norm([x-xopt,y-yopt,objfunc(x,y,FuncName)-objfunc(xopt,yopt,FuncName)]); % dICA(i+1) = d;
        dICA(i+1) = objfunc(xopt,yopt,FuncName) - objfunc(x,y,FuncName);
        % plot the current updated point
        FigCurrPoinPlot(nf,x0,y0,x,y,idx2d,i,OptAlg,d,FuncName)
        % reset (x0,y0) and switching index
        x0 = x; y0 = y; idx = 0;
    else
        % find the optimal step size following y coordinate
        alpopty = findoptalp(x0,y0,0,1,FuncName);
        % update value
        x = x0; 
        y = y0 + alpopty;
        % figure(nf + 2), hold on, plot(i,alpopty,'x');
        % compute the residue
        d = norm([x-xopt,y-yopt,objfunc(x,y,FuncName)-objfunc(xopt,yopt,FuncName)]); % dICA(i+1) = d;
        dICA(i+1) = objfunc(xopt,yopt,FuncName) - objfunc(x,y,FuncName);
        % plot the current updated point
        FigCurrPoinPlot(nf,x0,y0,x,y,idx2d,i,OptAlg,d,FuncName)
        % reset (x0,y0) and switching index
        x0 = x; y0 = y; idx = 1;
    end
end

% ==== optimization using new simultaneous coordinate ascent with line search ==== %

% define the maximization algorithm's name
OptAlg = 'SCM'; % (i.e., SCM algorithm)

% preset (x0,y0)
x0 = xinit; y0 = yinit;

% plot the objective function, initial point and optimal point
nf=5000; FigContMesh(nf,X,Y,F,x0,y0,xopt,yopt,idx2d,FuncName);

% preallocation of residues and first residue
dSCAn = zeros(1,niter+1); 
dSCAn(1) = norm([xinit-xopt,yinit-yopt,objfunc(xinit,yinit,FuncName)-objfunc(xopt,yopt,FuncName)]);

% run simultaneous coordinate ascent with line search and plot
for i = 1:niter
    % find the optimal step size with respect to x coordinate
    alpoptx = findoptalp(x0,y0,1,0,FuncName); 
    % find the optimal step size with respect to y coordinate
    alpopty = findoptalp(x0,y0,0,1,FuncName);
    % find the optimal step size for both x and y coordinates (eq. 5.24, page 97 of bouman's book)
    alpopt = findoptalp(x0,y0,alpoptx,alpopty,FuncName);
    % plot step-size
    figure(nf + 2), 
    subplot(2,1,1), hold on,
    if i<niter
        plot(i,alpoptx,'x','color',[0,0.5,0],'LineWidth',1); plot(i,alpopty,'o','color','b','LineWidth',1); 
    else
        plot(i,alpoptx,'x','color',[0,0.5,0],'LineWidth',1); plot(i,alpopty,'o','color','b','LineWidth',1);
        legend('a_{x}^{k}','a_{y}^{k}')
    end
    ylabel('a_{x}^{k}, a_{y}^{k}')
    subplot(2,1,2), hold on, 
    plot(i,alpopt,'*r'); ylabel('a_{0}^{k}')
    % update value
    x = x0 + alpopt*alpoptx;
    y = y0 + alpopt*alpopty;
    % compute the residue
    d = norm([x-xopt,y-yopt,objfunc(x,y,FuncName)-objfunc(xopt,yopt,FuncName)]); % dSCAn(i+1) = d;
    dSCAn(i+1) = objfunc(xopt,yopt,FuncName) - objfunc(x,y,FuncName);
    % plot the current updated point
    FigCurrPoinPlot(nf,x0,y0,x,y,idx2d,i,OptAlg,d,FuncName)
    % reset (x0,y0)
    x0 = x; y0 = y;
end

% ==== plot residues in function of iteration number ==== %
figure(6000)
semilogy(1:(niter+1),dGA,'-xr','LineWidth',1,'MarkerSize',7)
hold on
semilogy(1:(niter+1),dICA,'-o','LineWidth',1,'MarkerSize',6,'Color',[0,0.5,0],'MarkerFaceColor','w')
semilogy(1:(niter+1),dSCAn,'+-.','LineWidth',1,'MarkerSize',6,'Color',[0,0,0.5],'MarkerFaceColor','w')
legend({' GAFR',' ICMC',' SCM'},'FontSize',12)
legend('boxoff')
set(gca,'FontSize',12)
ylabel('convergence rate','FontSize',12), 
box off
