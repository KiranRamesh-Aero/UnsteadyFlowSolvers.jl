clc;clear all;close all
figure
set(gcf,'Units','Inches','Position',[1 1 3.25 2.6]);
set(gcf,'DefaultAxesFontName','Helvetica');
set(gcf,'DefaultTextFontName','Helvetica');
set(gcf,'DefaultAxesFontSize',8);
set(gcf,'DefaultTextFontSize',8);

set(gcf,'PaperUnits',get(gcf,'Units'));
pos = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 pos(3) pos(4)]);


theofull=load('force_mg_2d.dat');
%exp=load('exp_pitchramp_45_LE.dat');
%comp=load('comp_kd.dat');
comp=load('cl_mg_2d.csv');
exp=load('exp_mg_2d.csv');
cf=load('theo_mg_2d.csv');
%theofull_nosep=load('force_kd_nosep.dat');
%theofull_fullsep=load('force_kd_alwayssep.dat')
start_ind=floor(1*length(theofull(:,1))/2);
end_ind=length(theofull(:,1));
period_ind=floor(length(theofull(:,1))/2);
period=theofull(period_ind,1);

theo(:,:)=theofull(start_ind:end_ind,:);
%theo_nosep(:,:)=theofull_nosep(start_ind:end_ind,:);
%theo_fullsep(:,:)=theofull_fullsep(start_ind:end_ind,:);

theo(:,1)=(theo(:,1)-theofull(start_ind-1,1))/period;
%theo_nosep(:,1)=(theo_nosep(:,1)-theofull_nosep(start_ind-1,1))/period;
%theo_fullsep(:,1)=(theo_fullsep(:,1)-theofull_fullsep(start_ind-1,1))/period;

[ax, h1, h2] = plotyy(exp(:,1),exp(:,2),...
		      theo(:,1),theo(:,2));

set(ax(1),'YColor','k');
set(ax(2),'YColor','k');

set(h1,'Color','k');

set(h2,'Color',[0.5 0.5 0.5]);
set(h2,'LineWidth',1.2);

axes(ax(1));

hold on
plot(comp(:,1),comp(:,2),'-b','linewidth',0.1)
hcomp=plot(comp(1:18:end,1),comp(1:18:end,2),'.b','linewidth',1.5)
%plot(theo(:,1),smooth(medfilt1(theo(:,10),10)),'--','Color',[0 0.4 0],'linewidth',1.5) 
plot(theo(:,1),theo(:,10),'--','Color',[0 0.4 0],'linewidth',1.5)
plot(cf(:,1),cf(:,2),'r--','linewidth',1.5)
%plot(theo_nosep(:,1),smooth(theo_nosep(:,10)),'b--','linewidth',0.8)
%plot(theo_fullsep(:,1),smooth(smooth(theo_fullsep(:,10))),'r--','linewidth',0.8)

%hleg = legend('Experiment','CFD','LOM','Location','NorthWest');


%plot(0,interp1(theo(:,1),theo(:,10),0.01,'pchip'),'^k','markersize',5,'markerfacecolor','w','linewidth',1.5)
plot(0.1933,interp1(theo(:,1),theo(:,10),0.1695,'pchip'),'^k','markersize',4,'markerfacecolor','k','linewidth',1)
plot(0.2625,interp1(theo(:,1),theo(:,10),0.2625,'pchip'),'vk','markersize',4,'markerfacecolor','w','linewidth',1)
plot(0.5656,interp1(theo(:,1),theo(:,10),0.5656,'pchip'),'vk','markersize',4,'markerfacecolor','k','linewidth',1)
plot(0.7184,interp1(theo(:,1),theo(:,10),0.7184,'pchip'),'^k','markersize',4,'markerfacecolor','w','linewidth',1)
%plot(1.00,interp1(theo(:,1),theo(:,10),0.99,'pchip'),'^k','markersize',5,'markerfacecolor','k','linewidth',1.5)




%set(hleg,'Fontsize',7,'Box','off')
xlabel('t/T');
%ylabel('C_l')
  text(-0.1,5,'C_l');
set(gca,'XLim',[0 1],'YLim',[-30 30]);
set(gca,'XTick',[0:0.2:1],'YTick',[-30:10:30]);
% -grid on;

axes(ax(2));
set(gca,'YLim',[-60 60],'XLim',[0 1],'YTick',[-60:20:60]);
%h = text(4.4,41,'\alpha (right axis)');
%set(h,'Color',[0.5 0.5 0.5]);
%ylabel('\alpha, deg');
h = text(1.01,10,'\alpha (deg)');
set(h,'Color','k');

  

hold on;
plot(linspace(0.5,0.58,3),linspace(53,53,3),'.b','linewidth',1.5)
line([0.5 0.58],[53 53],'linewidth',0.1,'color','b')
     line([0.2 0.28],[53 53],'Color','k','linewidth',1.2)
text(0.29,53,'Exp','Fontsize',7)
    
line([0.2 0.28],[43 43],'Linestyle','--','Color',[0 0.4 0],'linewidth',1.5) 
text(0.29,43,'LDVM','Fontsize',7)
text(0.6,53,'CFL3D','Fontsize',7)
line([0.5 0.58],[43 43],'Linestyle','--','Color','r','linewidth',1.5) 
 text(0.6,43,'Theodorsen','Fontsize',7)
%line([0.03 0.11],[80 80],'Color',[0.5 0.5 0.5],'Linewidth',1.2)
%text(0.12,80,'\alpha','Fontsize',7)

text(0.72,28,'\alpha','Fontsize',8,'background','w')

%text(0.25,45,'LESP = \infty')
%text(0.52,-60,'LESP = 0')

set(ax(1),'Units','Inches');
set(ax(1),'Position',[0.35 0.45 2.5 2.0]);
set(ax(2),'Units','Inches');
set(ax(2),'Position',[0.35 0.45 2.5 2.0]);

print -depsc -loose ../figs/cl_mg_2d.eps

