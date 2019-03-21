clear all
close all
clc

global gsyn

init = zeros(4*2,1) + [.01 0 0 0 0 0 0 0].';
options = odeset('abstol',1e-8,'reltol',1e-8);
tspan = [0 30000];

gsyn = .020;
[t,Y] = ode45('USETHIS_TCdyn_multi',tspan,init,options);
g = length(Y)-2;
buf = Y(:,1);
while buf(g)<buf(g-1) || buf(g)<buf(g+1)
    g = g-1;
end
g1 = g;
g = g-200;
while buf(g)<buf(g-1) || buf(g)<buf(g+1)
    g = g-1;
end
g = g-200;
while buf(g)<buf(g-1) || buf(g)<buf(g+1)
    g = g-1;
end
g = g-200;
while buf(g)<buf(g-1) || buf(g)<buf(g+1)
    g = g-1;
end
g2 = g;
perorb_a = Y(g2:g1,:);
timeorb_a = t(g2:g1)-t(g2);


gsyn = .060;
[t,Y] = ode45('USETHIS_TCdyn_multi',tspan,init,options);
g = length(Y)-2;
buf = Y(:,1);
while buf(g)<buf(g-1) || buf(g)<buf(g+1)
    g = g-1;
end
g1 = g;
g = g-200;
while buf(g)<buf(g-1) || buf(g)<buf(g+1)
    g = g-1;
end
g = g-200;
while buf(g)<buf(g-1) || buf(g)<buf(g+1)
    g = g-1;
end
g = g-200;
while buf(g)<buf(g-1) || buf(g)<buf(g+1)
    g = g-1;
end
g2 = g;
perorb_s = Y(g2:g1,:);
timeorb_s = t(g2:g1)-t(g2);

subplot(2,2,1)
plot(timeorb_a,perorb_a(:,1)*100,'b','linewidth',2);
hold on
plot(timeorb_a,perorb_a(:,5)*100,'k','linewidth',2);
xlim([0 35])
ylabel('$V$ (mV)','interpreter','latex','fontsize',18)
set(gca,'xticklabel',[])
ylim([-80 0])
title('$\rho = 0.02\; {\rm mS}/{\rm cm}^2$','interpreter','latex','fontsize',16)

subplot(2,2,2)
theta = linspace(0,2*pi,2000);
x = cos(theta);
y = sin(theta);
plot(x,y,'--k');hold on
axis square
axis([-1.2 1.2 -1.2 1.2])
plot([0 1.2],[0 0],'k')
plot([0 1.6*cos(2)],[0 1.6*sin(2)],'k')
phasedif_a = 2*pi*7.898/15.8;
plot(cos(phasedif_a),sin(phasedif_a),'.b','markersize',30)
plot(cos(0),sin(0),'.k','markersize',30)




subplot(2,2,3)
plot(timeorb_s,perorb_s(:,1)*100,'b','linewidth',2);
hold on
plot(timeorb_s,perorb_s(:,5)*100,'k','linewidth',2);
phasedif_s = 2*pi*15.17/16.12;
ylabel('$V$ (mV)','interpreter','latex','fontsize',18)
xlabel('$t$ (ms)','interpreter','latex','fontsize',18)
ylim([-80 0])
xlim([0 35])
title('$\rho = 0.06 \;{\rm mS}/{\rm cm}^2$','interpreter','latex','fontsize',16)

subplot(2,2,4)
theta = linspace(0,2*pi,2000);
x = cos(theta);
y = sin(theta);
plot(x,y,'--k');hold on
axis square
axis([-1.2 1.2 -1.2 1.2])
plot([0 1.2],[0 0],'k')
plot([0 1.6*cos(2)],[0 1.6*sin(2)],'k')
plot(cos(phasedif_s),sin(phasedif_s),'.b','markersize',30)
plot(cos(0),sin(0),'.k','markersize',30)
text(0,0,'$\theta$','interpreter','latex','fontsize',16)
text(0,0,'$\theta$','interpreter','latex','fontsize',16)

text(0,0,'A','fontweight','bold','fontsize',20)
text(0,0,'B','fontweight','bold','fontsize',20)

addpath('C:\Users\dwilso81\Desktop\expfig')
% export_fig   -r600   -q101 modsim.png

return
