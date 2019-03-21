clear all
close all
clc

N = 1;
init = zeros(4*N,1)
options = odeset('abstol',1e-10,'reltol',1e-10);
tspan = [0 2000];

myfn = @TCdyn_single2;
[t,Y] = ode45(myfn,tspan,init,options);
plot(t,Y)

dt = .001;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
init = Y(length(Y),:);
tspan = [0:dt:30];
[t,Y] = ode45(myfn,tspan,init,options);


g = 2;buf = Y(:,1);
while buf(g)<buf(g-1) || buf(g)<buf(g+1)
    g = g+1;
end
g1 = g;
g = g+200;
while buf(g)<buf(g-1) || buf(g)<buf(g+1)
    g = g+1;
end
g2 = g;

perorb = Y(g1:g2,:);
timeorb = t(g1:g2)-t(g1);
init = perorb(1,:)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load inits2
% for k = 1:N
% sft = round(k*length(perorb)/N);
% init(1+(4*(k-1)):4+(4*(k-1))) = perorb(sft,1:4);
% end


T = timeorb(length(timeorb));
omega = 2*pi/T;
options = odeset('RelTol',3e-14,'AbsTol',3e-14);
Y = init;
Y(N*4+1) = T;Y = Y(:);
dy = ones(4,1);
while norm(dy)>.01
    
clear J
tspan = [0 Y(5)];
%%%Calculate Jacobian
for mm = 1:length(init)
eps = .01;
ptbplus = zeros(length(init),1);ptbplus(mm) = eps;
ptbminus = zeros(length(init),1);ptbplus(mm) = -eps;
cplus = Y(1:N*4) + ptbplus;
cminus = Y(1:N*4) +ptbminus;
[t,Yplus] = ode45(myfn,tspan,cplus,options);
[t,Yminus] = ode45(myfn,tspan,cminus,options);
mapplus = Yplus(length(Yplus),:);
mapminus = Yminus(length(Yminus),:);
J(:,mm) = [(mapplus-mapminus)/(2*eps)]';
end

J = J-eye(N*4);
%%%%Calculate Tdiff
epstime = 0.0001;
[t,Yplus] = ode45(myfn,[0 (Y(N*4+1)+epstime)],Y(1:N*4),options);
[t,Yminus] = ode45(myfn,[0 (Y(N*4+1)-epstime)],Y(1:N*4),options);
mapplus = Yplus(length(Yplus),:);
mapminus = Yminus(length(Yminus),:);
J(:,N*4+1) = [(mapplus-mapminus)/(2*epstime)]';
J(N*4+1,:) = [myfn(0,Y(1:4*N))' 0];

[t,Yorig] = ode45(myfn,[0 Y(N*4+1)],Y(1:N*4),options);
Yafter = Yorig(length(Yorig),:)';
b = [Y(1:N*4)-Yafter;0];

dy = J\b;
Y = Y+dy*.5;
[dy Y]
end

init  =Y(1:N*4)';tspan = 0:dt:Y(length(Y))
[timeorb,perorb] = ode45(myfn,tspan,init,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = timeorb(length(timeorb));
omega = 2*pi/T;

M = eye(length(init));
for k = 1:length(perorb)
    M = M + dt*numericaljacobian(myfn,perorb(k,:))*M;
end

[V,D] = eig(M);
VV = V;
DD = D;


[aa,bb] = sort(diag(D),'descend');
for k = 1:N-1
fa(k) = bb(k);
end



%%Get PRC
X1 = zeros(size(init));X1 = X1';
X1(1) = 1;
YY =  V.'
PRCinit = (YY\X1);
Zd = PRCinit;
myorb = fliplr(perorb.');
timeorb = flipud(timeorb);   

d = 1;
for mm = 1:1
for k = 1:length(myorb)-1
[k mm]
Zd(:,d+1) = numericaljacobian(myfn,myorb(:,k)).'*Zd(:,d)*dt - 0*eye(length(init))*Zd(:,d)*dt + Zd(:,d);
d = d+1;
end
end
Zd(:,1)./Zd(:,length(Zd))
mult = myfn(0,perorb(1,:))'*Zd(:,1)/omega;
Zd = Zd/mult;   %%%%NORMALIZE
PRC = Zd;
clear Zd;







PRC = fliplr(PRC);
timeorb = flipud(timeorb);

PRC = real(PRC)


% save('PRC_synaptic_gsynp015_esyn_n56','perorb','PRC','timeorb','T')









return