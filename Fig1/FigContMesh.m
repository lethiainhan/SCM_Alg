function FigContMesh(nf,X,Y,F,x0,y0,xopt,yopt,idx2d,FuncName)
% plot the contour or mesh of objective function, initial point and optimal point
% nf is the figure number 1000, 2000, 3000, 4000, ...
% idx2d is index for 2D or 3D plot: 1 = 2d and 0 = 3d

if idx2d == 1
    figure(nf); contour(X,Y,F,50,'HandleVisibility','off'); hold on; 
    plot(x0,y0,'o','LineWidth',1,'MarkerSize',6,'MarkerEdgeColor',[0,0.5,0]);
    plot(xopt,yopt,'+','LineWidth',1,'MarkerSize',7,'MarkerEdgeColor',[1,0.5,0]);
    legend({' initiation',' optimum'},'FontSize',12,'Location','southwest')
    legend('boxoff')
    xlabel('x','FontSize',12), ylabel('y','FontSize',12)
    set(gca,'FontSize',12)
else
    figure(nf+1); mesh(X,Y,F); view(30,20); colorbar; hold on; 
    plot3(x0,y0,objfunc(x0,y0,FuncName),'s','LineWidth',1,'MarkerSize',8,'MarkerEdgeColor',[0,0.5,0]);
    plot3(xopt,yopt,objfunc(xopt,yopt,FuncName),'o','LineWidth',1,'MarkerSize',7,'MarkerEdgeColor',[1,0.5,0]);
    legend({' obj. func.',' init. val.',' opti. val.'},'FontSize',12,'Location','northwest')
    xlabel('x','FontSize',12), ylabel('y','FontSize',12)
    set(gca,'FontSize',12)
end

end