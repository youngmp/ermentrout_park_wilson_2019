clear all
close all
clc

global F1 F2 F3 F4 F5 F6 theta alpha T1 T2 F1_t1mt2 F2_t1mt2 F3_t1mt2 F4_t1mt2 F5_t1mt2 F6_t1mt2 kappas

load reducedparameters_noselfcoupling_dtp0001_Ib3p75

omega = 2*pi/T
V = perorb(:,1);
S = perorb(:,4);
Z = PRC(1,:);Z = Z(:);   %Z = Z/omega;
B = B{1}(1,:);B = B(:);  %B = B/omega;
C = C{1}(1,:);C = C(:);
I = IRC{1}(1,:);I = I(:);
pv = ydiffall{1}(1,:);pv = pv(:);
ps = ydiffall{1}(4,:);ps = ps(:);
theta = linspace(0,2*pi,length(V));theta = theta(:);
T = max(timeorb);
pr = ydiffall{1}(3,:);pr = pr(:);

subplot(2,2,1);hold on
plot(theta,Z)
subplot(2,2,2);hold on
plot(theta,I)
subplot(2,2,3);hold on
plot(theta,B)
subplot(2,2,4);hold on
plot(theta,C)




theta2 = linspace(0,2*pi,1000);
V = interp1(theta,V,theta2).';
S = interp1(theta,S,theta2).';
Z = interp1(theta,Z,theta2).';
B = interp1(theta,B,theta2).';
C = interp1(theta,C,theta2).';
I = interp1(theta,I,theta2).';
pv = interp1(theta,pv,theta2).';
ps = interp1(theta,ps,theta2).';
theta = theta2;


alpha = -.015;  %%i_syn comes in with a minus sign
esyn = -60;


for k = 1:length(Z)
        f1 = circshift(Z,-k);
        f3 = circshift(V,-k);
        f5 = circshift(S,-k);
        f6 = circshift(B,-k);
        f7 = circshift(pv,-k);
        f9 = circshift(C,-k);
        f10 = circshift(I,-k);
        f11 = circshift(ps,-k);
        
        %%phase equation
        F1_t1mt2(k) = mean(f1.*(S + f5).*(f3-esyn));
        F2_t1mt2(k) = mean(f6.*(f5+S).*(f3-esyn) + f1.*f11.*(f3-esyn) + f1.*(f5+S).*f7);
        F3_t1mt2(k) = mean(f1.*ps.*(f3-esyn));
        
        %%isostable equations
        F4_t1mt2(k) = mean(f10.*(S + f5).*(f3-esyn));
        F5_t1mt2(k) = mean(f9.*(f5+S).*(f3-esyn) + f10.*(f11).*(f3-esyn)  + f10.*(f5+S).*f7 );
        F6_t1mt2(k) = mean( f10.*ps.*(f3-esyn) );
end
%  F1_t1mt2 = fliplr(F1_t1mt2);
%  F2_t1mt2 = fliplr(F2_t1mt2);
%  F3_t1mt2 = fliplr(F3_t1mt2);
%  F4_t1mt2 = fliplr(F4_t1mt2);
%  F5_t1mt2 = fliplr(F5_t1mt2);
%  F6_t1mt2 = fliplr(F6_t1mt2);

F1_t1mt2 = [F1_t1mt2 F1_t1mt2];
F2_t1mt2 = [F2_t1mt2 F2_t1mt2];
F3_t1mt2 = [F3_t1mt2 F3_t1mt2];
F4_t1mt2 = [F4_t1mt2 F4_t1mt2];
F5_t1mt2 = [F5_t1mt2 F5_t1mt2];
F6_t1mt2 = [F6_t1mt2 F6_t1mt2];


% FF = .5*mean(F1_t1mt2)*ones(size(theta));
% for k = 1:100
%     b(k) = mean(sin(k*pi*theta/(2*pi)).*F1_t1mt2);
%     a(k) = mean(cos(k*pi*theta/(2*pi)).*F1_t1mt2);
%     FF = FF + a(k)*cos(k*pi*theta/(2*pi)) + b(k)*sin(k*pi*theta/(2*pi))
% end
% plot(FF);hold on;plot(F1_t1mt2)


%alpha should have a minus sign in this derivation
alpha = -.06;
init = [0 0 0]';
for k = 1:10
J = numericaljacobian(@isostablereducedmodel_better_firstorder,init);
buf = isostablereducedmodel_better_firstorder(0,init);
dx = -J\buf;
init = init + dx;
end
eig(J)
init


return

figure
subplot(1,2,1)
buf = -(F1_t1mt2 - fliplr(F1_t1mt2));
plot(linspace(-2*pi,2*pi,length(buf)),buf,'k','linewidth',2);hold on
plot([0 2*pi],[0 0],'--k')
xlim([0 2*pi])
subplot(1,2,2)
plot(linspace(-2*pi,2*pi,length(buf)),buf,'k','linewidth',2);hold on
plot([-pi pi],[0 0],'--k')
xlim([-.65 .65])
d = 1;
for alpha = -.01:-.03:-.2
init = [0 0 0]';
%alpha should have a minus sign in this derivation
for k = 1:20
J = numericaljacobian(@isostablereducedmodel_better,init);
buf = isostablereducedmodel_better(0,init);
dx = -J\buf;
init = init + dx;
end
initov(:,d) = init;d = d+1;
hold on
buf = -(F1_t1mt2 - fliplr(F1_t1mt2)  +   init(2)*(F2_t1mt2 - fliplr(F3_t1mt2)) +  init(3)*(F3_t1mt2 - fliplr(F2_t1mt2)));
subplot(1,2,1); hold on
plot(linspace(-2*pi,2*pi,length(buf)),buf);
subplot(1,2,2); hold on
plot(linspace(-2*pi,2*pi,length(buf)),buf);
end

subplot(1,2,1)
ylabel('$\dot{\Upsilon} / \rho$','interpreter','latex','fontsize',18)
xlabel('$\Upsilon$','interpreter','latex','fontsize',18)
ylim([-.3 .3])

subplot(1,2,1)
xlabel('$\Upsilon$','interpreter','latex','fontsize',18)
ylim([-.025 .025])


text(0,0,'A','fontweight','bold','fontsize',20)
text(0,0,'B','fontweight','bold','fontsize',20)
text(0,0,'C','fontweight','bold','fontsize',20)

text(0,0,'increasing $\rho$','interpreter','latex','fontsize',18)
% text(0,0,'$\rho$','interpreter','latex','fontsize',18)

subplot(5,5,25)
box on
plot(linspace(-2*pi,2*pi,length(buf)),-F4_t1mt2,'k','linewidth',2);hold on
ylabel('$H_4(\Upsilon)$','interpreter','latex','fontsize',18)
xlim([0 2*pi])

set(gca,'xtick',[0 pi 2*pi])
xt1{1} = 0;
xt1{2} = '\pi';
xt1{3} = '2\pi';
set(gca,'xticklabel',xt1)

addpath('C:\Users\dwilso81\Desktop\expfig')
% export_fig   -r600   -q101 bifnfig.png




return


figure
alpha = -.06
options = odeset('abstol',1e-10,'reltol',1e-10,'maxstep',3);
init = init + [+.0001 0 0]';
[t,Y] = ode45(@isostablereducedmodel_better,[0 30000],init,options)
clf
subplot(1,2,2)
plot(t,Y(:,1))
hold on
subplot(1,2,1)
plot(t,Y(:,2));
hold on
plot(t,Y(:,3))
init = Y(length(Y),:).';
init2 = init;

amin = -.04
amax = -.06;
for j = 1:20
alpha = (amin+amax)/2

for k = 1:40
J = numericaljacobian(@isostablereducedmodel_better,init);
buf = isostablereducedmodel_better(0,init);
dx = -J\buf;
init = init + dx;
end

if norm(init2-init)<1  &&  norm(dx)<.001
    init2 = init;
    amax = alpha;
else
    amin = alpha;
end

end



init = init2;
alpha = amax;
for k = 1:40
J = numericaljacobian(@isostablereducedmodel_better,init);
buf = isostablereducedmodel_better(0,init);
dx = -J\buf;
init = init + dx;
end

subplot(1,2,1);
buf = -(F1_t1mt2 - fliplr(F1_t1mt2));
plot(linspace(-2*pi,2*pi,length(buf)),buf,'k','linewidth',2);
xlim([-pi pi])
buf = -(F1_t1mt2 - fliplr(F1_t1mt2)  +   init(2)*(F2_t1mt2 - fliplr(F3_t1mt2)) +  init(3)*(F3_t1mt2 - fliplr(F2_t1mt2)));
subplot(1,2,1); hold on
plot(linspace(-2*pi,2*pi,length(buf)),buf);
subplot(1,2,2); hold on
plot(linspace(-2*pi,2*pi,length(buf)),buf);









return
