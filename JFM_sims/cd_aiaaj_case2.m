
figure
set(gcf,'Units','Inches','Position',[1 1 3.25 2.6]);
set(gcf,'DefaultAxesFontName','Helvetica');
set(gcf,'DefaultTextFontName','Helvetica');
set(gcf,'DefaultAxesFontSize',8);
set(gcf,'DefaultTextFontSize',8);

set(gcf,'PaperUnits',get(gcf,'Units'));
pos = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 pos(3) pos(4)]);


theo=load('force_cfdjfm_2.dat');
%exp=load('exp_pitchramp_45_LE.dat');
comp=load('comp_cfdjfm_2.dat');

[ax, h1, h2] = plotyy(comp(:,1),comp(:,6),...
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
plot(theo(:,1),theo(:,11),'--','Color',[0 0.4 0],'linewidth',1.5) 

plot(1.0160,interp1(theo(:,1),theo(:,11),1.0160,'pchip'),'vk','markersize',4,'markerfacecolor','w','linewidth',1)
plot(1.2418,interp1(theo(:,1),theo(:,11),1.2418,'pchip'),'vk','markersize',4,'markerfacecolor','k','linewidth',1)
plot(6.1285,interp1(theo(:,1),theo(:,11),6.1285,'pchip'),'vk','markersize',4,'markerfacecolor','w','linewidth',1)
%plot(8.0,interp1(theo(:,1),theo(:,11),8.0,'pchip'),'vk','markersize',4,'markerfacecolor','k','linewidth',1)
plot(2.1450,interp1(theo(:,1),theo(:,11),2.1450,'pchip'),'^k','markersize',4,'markerfacecolor','w','linewidth',1)
plot(5.3544,interp1(theo(:,1),theo(:,11),5.3544,'pchip'),'^k','markersize',4,'markerfacecolor','k','linewidth',1)


%hleg = legend('Experiment','CFD','LOM','Location','NorthWest');

%set(hleg,'Fontsize',7,'Box','off')
xlabel('t^*');
%ylabel('C_l')
  text(-0.7,3.5,'C_d');
set(gca,'XLim',[0 8],'YLim',[-1 5]);
set(gca,'XTick',[0:1:8],'YTick',[-1:1:5]);
% -grid on;

axes(ax(2));
set(gca,'YLim',[0 90],'XLim',[0 8],'YTick',[0:15:90]);
%h = text(4.4,41,'\alpha (right axis)');
%set(h,'Color',[0.5 0.5 0.5]);
%ylabel('\alpha, deg');
h = text(8.05,67.5,'\alpha (deg)');
set(h,'Color','k');

  

hold on;
%plot(linspace(1.9,2.9,8),linspace(56,56,8),'.b','linewidth',1.5)
%line([1.9 2.9],[56 56],'linewidth',0.1,'color','b')
     line([0.1 0.8],[85 85],'Color','k')
text(0.9,85,'CFD','Fontsize',7)
    
%text(3.0,56,'CFD','Fontsize',7)
 line([5.1 5.8],[85 85],'Linestyle','--','Color',[0 0.4 0],'linewidth',1.5) 
text(5.9,85,'LDVM','Fontsize',7)
%line([0.1 0.8],[85 85],'Color',[0.5 0.5 0.5],'Linewidth',1.2)
%text(0.9,85,'\alpha','Fontsize',7)
%text(1.5,30,'\alpha','Fontsize',8,'background','w')
%text(2.5,75,'\alpha','Fontsize',8,'background','w')
%text(2.4,70,'\alpha','Fontsize',8,'background','w')
%set(gcf,'color','white')
%set(gcf,'inverthardcopy','off')
text(2.8,81,'\alpha','Fontsize',8,'background','w')


set(ax(1),'Units','Inches');
set(ax(1),'Position',[0.35 0.45 2.5 2.0]);
set(ax(2),'Units','Inches');
set(ax(2),'Position',[0.35 0.45 2.5 2.0]);

print -depsc -loose ../figs/cd_aiaaj_case2.eps

