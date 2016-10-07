
figure
set(gcf,'Units','Inches','Position',[1 1 3.25 2.6]);
set(gcf,'DefaultAxesFontName','Helvetica');
set(gcf,'DefaultTextFontName','Helvetica');
set(gcf,'DefaultAxesFontSize',8);
set(gcf,'DefaultTextFontSize',8);

set(gcf,'PaperUnits',get(gcf,'Units'));
pos = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 pos(3) pos(4)]);


theo=load('force_eld_90_lesp0.11.dat');
%exp=load('exp_pitchramp_45_LE.dat');
comp=load('comp_eld_90.dat');

[ax, h1, h2] = plotyy(comp(:,1),comp(:,12)*2,...
		      theo(:,1),theo(:,2));

set(ax(1),'YColor','k');
set(ax(2),'YColor','k');

set(h1,'Color','k');

set(h2,'Color',[0.5 0.5 0.5]);
set(h2,'LineWidth',1.2);

axes(ax(1));

hold on
%plot(comp(:,1),comp(:,2),'-b','linewidth',0.1)
%hcomp=plot(comp(1:18:end,1),comp(1:18:end,2),'.b','linewidth',1.5)
plot(theo(:,1),smooth(theo(:,10)),'--','Color',[0 0.4 0],'linewidth',1.5) 

plot(1.2110,interp1(theo(:,1),theo(:,10),1.211,'pchip'),'^k','markersize',4,'markerfacecolor','w','linewidth',1)
%plot(5,interp1(theo(:,1),theo(:,10),5,'pchip'),'^k','markersize',5,'markerfacecolor','k','linewidth',1.5)


%hleg = legend('Experiment','CFD','LOM','Location','NorthWest');

%set(hleg,'Fontsize',7,'Box','off')
xlabel('t^*');
%ylabel('C_l')
  text(-0.5,2.5,'C_l');
set(gca,'XLim',[0 5],'YLim',[-1 5]);
set(gca,'XTick',[0:1:5],'YTick',[-1:1:5]);
% -grid on;

axes(ax(2));
set(gca,'YLim',[0 90],'XLim',[0 5],'YTick',[0:15:90]);
%h = text(4.4,41,'\alpha (right axis)');
%set(h,'Color',[0.5 0.5 0.5]);
%ylabel('\alpha, deg');
h = text(5.05,52.5,'\alpha (deg)');
set(h,'Color','k');

  

hold on;
%plot(linspace(1.9,2.9,8),linspace(56,56,8),'.b','linewidth',1.5)

    line([0.1 0.8],[85 85],'Color','k')
text(0.9,85,'CFD','Fontsize',7)
    
%text(3.0,56,'CFD','Fontsize',7)
 line([1.7 2.4],[85 85],'Linestyle','--','Color',[0 0.4 0],'linewidth',1.5) 
text(2.5,85,'LDVM','Fontsize',7)
%line([0.1 0.8],[85 85],'Color',[0.5 0.5 0.5],'Linewidth',1.2)
%text(0.9,85,'\alpha','Fontsize',7)

%text(2.5,37,'\alpha','Fontsize',8,'background','w')
text(2.5,37,'\alpha','Fontsize',8,'background','w')



set(ax(1),'Units','Inches');
set(ax(1),'Position',[0.35 0.45 2.5 2.0]);
set(ax(2),'Units','Inches');
set(ax(2),'Position',[0.35 0.45 2.5 2.0]);

print -depsc -loose ../figs/cl_eld_90.eps

