clear all
close all
clc

global F1 F2 F3 F4 F5 F6 theta alpha T1 T2 F1_t1mt2 F2_t1mt2 F3_t1mt2 F4_t1mt2 F5_t1mt2 F6_t1mt2 kappas


% theta1 = mod(y(1),2*pi);
% theta2 = mod(y(2),2*pi);
% 
% dy = zeros(2,1);
% 
% 
% dy(1) = qinterp2(T1,T2,F1,theta2,theta1);
% dy(2) = qinterp2(T1,T2,F1,theta1,theta2);

load reducedparameters_noselfcoupling

V = perorb(:,1);
S = perorb(:,4);
Z = PRC(1,:);Z = Z(:);
Esyn = -61;
B = B{1}(1,:);B = B(:);
C = C{1}(1,:);C = C(:);
I = IRC{1}(1,:);I = I(:);
pv = ydiffall{1}(1,:);pv = pv(:);
ps = ydiffall{1}(4,:);ps = ps(:);
theta = linspace(0,2*pi,length(V));theta = theta(:);
T = max(timeorb);
pr = ydiffall{1}(3,:);pr = pr(:);



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
esyn = -61;

F1 = zeros(1000,1000);
F2 = zeros(1000,1000);
F3 = zeros(1000,1000);
F4 = zeros(1000,1000);
F5 = zeros(1000,1000);
F6 = zeros(1000,1000);

for k = 1:length(Z)
    for j = 1:length(Z)
        f1 = circshift(Z,-k);
        f2 = circshift(V,-j);
        f3 = circshift(V,-k);
        f4 = circshift(S,-j);
        f5 = circshift(S,-k);
        
        f6 = circshift(B,-k);
        f7 = circshift(pv,-k);
        f8 = circshift(ps,-j);
        f9 = circshift(C,-k);
        f10 = circshift(I,-k);
        f11 = circshift(ps,-k);
        
        %%phase equation
        F1(k,j) = mean(f1.*(f4 + f5).*(f3-esyn));
        F2(k,j) = mean(f6.*f4.*(f3-esyn) + f1.*(f4.*f7) + f1.*f11.*(f3-esyn));
        F3(k,j) = mean(f1.*f8.*(f3-esyn));
        
        %%isostable equations
        F4(k,j) = mean(f10.*(f4 + f5).*(f3-esyn));
        F5(k,j) = mean(f9.*(f4).*(f3-esyn) + f10.*(f4.*f7)  + f10.*f11.*(f3-esyn)  );
        F6(k,j) = mean( f10.*f8.*(f3-esyn) );
    end
    k
    
    F1_t1mt2(k) = F1(k,j);   %%maybe these need to be flipped...
    F2_t1mt2(k) = F2(k,j);
    F3_t1mt2(k) = F3(k,j);
    F4_t1mt2(k) = F4(k,j);
    F5_t1mt2(k) = F5(k,j);
    F6_t1mt2(k) = F6(k,j);
    
end

 F1_t1mt2 = fliplr(F1_t1mt2);
 F2_t1mt2 = fliplr(F2_t1mt2);
 F3_t1mt2 = fliplr(F3_t1mt2);
 F4_t1mt2 = fliplr(F4_t1mt2);
 F5_t1mt2 = fliplr(F5_t1mt2);
 F6_t1mt2 = fliplr(F6_t1mt2);

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




[T1,T2] = meshgrid(theta,theta);
surf(T1,T2,F1,'edgecolor','none');



%alpha should have a minus sign in this derivation
alpha = -.01;
init = [0 0 0]';
for k = 1:10
J = numericaljacobian(@isostablereducedmodel_better_firstorder,init);
buf = isostablereducedmodel_better_firstorder(0,init);
dx = -J\buf
init = init + dx;
end
eig(J)

return
init = [0 0 0]';
% init = Y(length(Y),:).';
%alpha should have a minus sign in this derivation
alpha = -.001;
for k = 1:20
J = numericaljacobian(@isostablereducedmodel_better,init);
buf = isostablereducedmodel_better(0,init);
dx = -J\buf;
init = init + dx;
end
eig(J)
init


clf
subplot(1,2,2)
buf = alpha*(F1_t1mt2 - fliplr(F1_t1mt2)  +   init(2)*(F2_t1mt2 - fliplr(F3_t1mt2)) +   init(3)*(F3_t1mt2 - fliplr(F2_t1mt2)));
plot(linspace(-2*pi,2*pi,length(buf)),buf);
hold on
buf = alpha*(F1_t1mt2 - fliplr(F1_t1mt2));
plot(linspace(-2*pi,2*pi,length(buf)),buf)
subplot(1,2,1)
plot(linspace(-2*pi,2*pi,length(buf)),buf)


options = odeset('abstol',1e-10,'reltol',1e-10,'maxstep',3);
init = init + [.0001 0 0]';
[t,Y] = ode45(@isostablereducedmodel_better,[0 2000],init,options)

clf
subplot(1,2,2)
plot(t,Y(:,1))
hold on
subplot(1,2,1)
plot(t,Y(:,2));
hold on
plot(t,Y(:,3))



return
