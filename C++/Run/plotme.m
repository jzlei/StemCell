function plotme()
global A B
path='output/';
A=load(strcat(path,'md.dat'));

fig=figure(1);
clf();

%subplot(2,2,1);
hold on;
%line([1666.7 1666.7],[0 7000],'linestyle','--','color','k');
plot((A(1:240:end,1)-10000)/24,A(1:240:end,3)/1000);

xlabel('Day');
ylabel('Cells (x 10^3)');
box on;
set(gca,'fontsize',12);
xlim([0 3500]);
ylim([4200 5600]/1000);

T=[2000,20000,45000,60000,80000];

x=[0.18,0.35,0.60,0.55,0.80];
y=[0.54,0.60,0.20,0.79,0.42];

y0=A((10000+T)*4,3)/1000;

x1=[350,1100,2100,2200,3200];
y1=[4880,4990,4500,5320,4920]/1000;

for i=1:5
    plot(T(i)/24,y0(i),'ko','markerfacecolor','k','markersize',5);
    plot([T(i)/24,x1(i)],[y0(i),y1(i)],'k-','linewidth',0.5);
end

for i=1:5
    B=load(strcat(path,'md-',num2str((10000+T(i))*4),'.dat'));
    axes('position',[x(i),y(i),0.1,0.1]);
    hold on;
    c=findcolor(B);
    plot([0 1],[0.3 0.3],'k--','linewidth',0.5);
    plot([0.3 0.3],[0 1],'k--','linewidth',0.5);
%    plot(B(:,2),B(:,3),'b.','markersize',1);
    scatter(B(:,2),B(:,3),1,c,'filled');
    xlabel('x1');
    ylabel('x2');
    set(gca,'xticklabel',{''});
    set(gca,'yticklabel',{''});
    box off;
    xlim([0 1]);
    ylim([0 1]);
    set(gca,'fontsize',8);
end

exportfig(fig,'fig8.eps','color','cmyk','fontmode','scaled','fontsize',1);

end

function Z=density(A)
dx=0.01;
Z=zeros(101,101);
n=size(A,1);
for k=1:n
    i=floor(A(k,2)/dx)+1;
    j=floor(A(k,3)/dx)+1;
    Z(i,j)=Z(i,j)+1;
end
end

function c=findcolor(A)
dx=0.01;
n=size(A,1);
Z=zeros(101,101);
c=zeros(n,1);
for k=1:n
    i=floor(A(k,2)/dx)+1;
    j=floor(A(k,3)/dx)+1;
    Z(i,j)=Z(i,j)+1;
end
for k=1:n
    i=floor(A(k,2)/dx)+1;
    j=floor(A(k,3)/dx)+1;
    c(k)=Z(i,j)/n;
end
end