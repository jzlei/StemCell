function plottissuegrowth()
global A B
A=load('output/md-ntpx.dat');
B=load('output/md.dat');

m=70;
t0=B(m,1)/24;
col=[0.2,0.2,1;
    0.2,1,0.2;
    1,0.2,0.2;
    0.1,0.1,0.1];
S=A(:,4)+A(:,5);

fig=figure(1);
clf();
subplot(4,2,1);
loglog(A(:,1)/24,S,'k');
hold on;
%plot([t0 t0],[1e-4 1e8],'--','color',[0.7,0.7,0.7]);
xlim([1e0 1e4]);
ylim([5e1 1e6]);
set(gca,'ytick',[1e2,1e3,1e4,1e5]);
set(gca,'xticklabel',{''});
ylabel('Cells');
text(2e-1,1e6,'A','fontweight','bold');


axes('position',[0.13 0.585 0.335 0.165]);
semilogx(A(:,1)/24,100*A(:,2)./S,'color',col(1,:));
hold on;
semilogx(A(:,1)/24,100*A(:,3)./S,'color',col(2,:));
semilogx(A(:,1)/24,100*A(:,5)./S,'color',col(3,:));
%plot([t0 t0],[-5 105],'--','color',[0.7,0.7,0.7]);
xlim([1e0 1e4]);
ylim([-5 105]);
set(gca,'xtick',[1e0 1e1 1e2 1e3 1e4]);
xlabel('Day');
ylabel('Cells (%)');
text(5e2,65,'Stem cells','color',col(2,:),'fontsize',8);
text(5e2,53,'Progenitor cells','color',col(1,:),'fontsize',8);
text(5e2,41,'TD cells','color',col(3,:),'fontsize',8);

subplot(2,2,2);
hold on;
[m,n]=size(B);
m=70;
x=0:0.01:1;
Y=B(m-1,2:n)+B(m,2:n);
Y0=sum(Y);
Z=Y/Y0;
plot(x,B(m-1,2:n)/Y0,'color',(col(1,:)+col(2,:))/2);
plot(x,B(m,2:n)/Y0,'color',col(3,:));
plot(x,Z, 'color',col(4,:));
legend('SP cells','TD cells','SP + TD cells');
xlim([0 1]);
ylim([0 0.04]);
box on;
xlabel('x');
ylabel('Density');
text(-0.22,0.03, 'B', 'fontweight','bold');
exportfig(fig,'figures/tissue.eps','color','cmyk','fontmode','scaled','fontsize',1);
end