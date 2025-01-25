% This script is used to plot the 2D trid function and to evaluate the average 
% computational time to reach an epsi-optimal solution of with respect to the 
% number of variables when three SCM, GAFR and CICM methods are respectively used

% ------------------------ Plot 2D trid function ------------------------ %

% define the rectangular grid [-J^2,J^2]x[-J^2,J^2] = [-4,4]x[-4,4]
m = 101; [X,Y] = meshgrid(linspace(-4,4,m),linspace(-4,4,m));

% compute mesh of objectif junction
F = -(X-1).^2 - (Y-1).^2 - X.*Y;

% mesh plot
figure; meshc(X,Y,F); 
xlabel('x','FontSize',16), ylabel('y','FontSize',16)
zlabel('F_2(x,y)','FontSize',16)
set(gca,'FontSize',16)
view(25,16)

% ----------------- Evaluate average computational time ----------------- %

% set profile on to turn off JIT compilation in matlab for "fair" comparison
% (see e.g., https://lips.cs.princeton.edu/jit-compilation-in-matlab/)
clear; profile off; profile clear; profile on;

% desired acceptable error to the maximum
epsi = 1e-6;

% number of repetitions to compute the average quantities
Nr = 100;

% vector of number of variable (i.e., dimension of the trid function)
vJ = [10:10:90, 100:50:1000];

% pre-allocation for vectors of average computation time
nJ = numel(vJ); vETscm = zeros(1,nJ); vETgafr = vETscm; vETcicm = vETscm;

% main loop for computation
for n = 1:nJ
    
    % display the current loop 
    fprintf('iteration n°: %d \n',n);
    
    % compute the average computation time
    [vETscm(n),vETgafr(n),vETcicm(n)] = ExecTimeMaxAlgTridFunc(vJ(n),epsi,Nr);
end

% plot the average computational time w.r.t. the number of variables
figure;
hold on
plot(vJ,vETscm,'-ro'); % SCM method
plot(vJ,vETgafr,'-bx'); % GAFR method
plot(vJ,vETcicm,'-k*'); % CICM method
legend({' SCM',' GAFR',' CICM'},'FontSize',12,'Location','northwest')
legend('boxoff')
set(gca,'FontSize',12)
xlabel('n° of variables','FontSize',12), 
ylabel('average computational time (s)','FontSize',12), 
box off

% set profile off to turn on JIT compilation in matlab for normal use
profile off;
