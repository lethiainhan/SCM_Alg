function FigCurrPoinPlot(nf,x0,y0,x,y,idx2d,i,OptAlg,d,FuncName)
% plot the current point in the contour or mesh of objective function
% nf is the figure number 1000, 2000, 3000, 4000, ...
% idx2d is index for 2D or 3D plot: 1 = 2d and 0 = 3d
% OptAlg is the name of algorithm among 'GA', 'ICA', 'SCA_o','SCA_n'

if idx2d == 1
    figure(nf); 
    quiver(x0,y0,x-x0,y-y0,0,'MaxHeadSize',0.75,'HandleVisibility','off','color','r','LineWidth',0.75)
else
    figure(nf+1); plot3([x0 x],[y0, y],[objfunc(x0,y0,FuncName), objfunc(x,y,FuncName)],'-*r','LineWidth',0.75,'MarkerSize',6,'HandleVisibility','off')
    title([OptAlg,' - n° iter = ', num2str(i),', residue = ',num2str(d)],'FontSize',12,'FontWeight','normal')
end

end