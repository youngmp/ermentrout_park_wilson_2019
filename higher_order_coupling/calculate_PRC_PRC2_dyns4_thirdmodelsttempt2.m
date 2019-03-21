clear all
close all
clc

N = 1;
init = zeros(4*N,1)
options = odeset('abstol',1e-12,'reltol',1e-12);
tspan = [0 2000];

myfn = @TCdyn_single2;
[t,Y] = ode45(myfn,tspan,init,options);
plot(t,Y)

dt = .00005;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
init = Y(length(Y),:);
tspan = [0:dt:50];
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
perorb = circshift(perorb,round(length(perorb)/4))

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
Y(N*length(init)+1) = T;Y = Y(:);
dy = size(init);
while norm(dy)>.001
    
clear J
tspan = [0 Y(5)];
%%%Calculate Jacobian
for mm = 1:length(init)
eps = .01;
ptbplus = zeros(length(init),1);ptbplus(mm) = eps;
ptbminus = zeros(length(init),1);ptbplus(mm) = -eps;
cplus = Y(1:N*length(init)) + ptbplus;
cminus = Y(1:N*length(init)) +ptbminus;
[t,Yplus] = ode45(myfn,tspan,cplus,options);
[t,Yminus] = ode45(myfn,tspan,cminus,options);
mapplus = Yplus(length(Yplus),:);
mapminus = Yminus(length(Yminus),:);
J(:,mm) = [(mapplus-mapminus)/(2*eps)]';
end

J = J-eye(N*length(init));
%%%%Calculate Tdiff
epstime = 0.0001;
[t,Yplus] = ode45(myfn,[0 (Y(N*length(init)+1)+epstime)],Y(1:N*length(init)),options);
[t,Yminus] = ode45(myfn,[0 (Y(N*length(init)+1)-epstime)],Y(1:N*length(init)),options);
mapplus = Yplus(length(Yplus),:);
mapminus = Yminus(length(Yminus),:);
J(:,N*length(init)+1) = [(mapplus-mapminus)/(2*epstime)]';
J(N*length(init)+1,:) = [myfn(0,Y(1:length(init)*N))' 0];

[t,Yorig] = ode45(myfn,[0 Y(N*length(init)+1)],Y(1:N*length(init)),options);
Yafter = Yorig(length(Yorig),:)';
b = [Y(1:N*length(init))-Yafter;0];

dy = J\b;
Y = Y+dy;
[dy Y]
end

init  =Y(1:N*length(init))';tspan = 0:dt:Y(length(Y))
[timeorb,perorb] = ode45(myfn,tspan,init,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = timeorb(length(timeorb));
omega = 2*pi/T;

init = perorb(1,:);
M = eye(length(init));
for k = 1:length(perorb)
  test1 = numericaljacobian(@TCdyn_single2,perorb(k,:));
  M = M + dt*(test1)*M;
end

[VV,DD ] = eig(M);
D = DD;
V = VV;

[aa,bb] = sort(diag(D),'descend')
fa = bb(2)
lambda = (D(fa,fa)); lambda = diag(lambda)
kappa = log(lambda)/T;
omega = 2*pi/T;






%%%%%%Calculate ISOSTABLE Response Curves
%%Calculate P
psinot = 0.0001;
for k = 1:length(fa)
init2 = init + V(:,fa(k)).'*psinot;
[tu,Yu] = ode113('TCdyn_single2',[0:dt:T],init,options);
perorb = Yu;
timeorb = tu;
timeorb = flipud(timeorb);
[tp,Yp] = ode113('TCdyn_single2',[0:dt:T],init2,options);
ydiff{k} = (Yp-Yu)/psinot;
((ydiff{k}(length(ydiff{k}),:)))./(ydiff{k}(1,:))
end



[dum,www] = min(abs(diag(DD)-1))

%%Get PRC
X1 = [0 0 0 0].';
X1(www) = 2;
YY =  V.'
IRCinit = (YY\X1);
Zd = IRCinit;
myorb = fliplr(perorb.');

d = 1;
for k = 1:length(myorb)-1
[k mm]
Zd(:,d+1) = numericaljacobian(@TCdyn_single2,myorb(:,k)).'*Zd(:,d)*dt - 0*eye(length(init))*Zd(:,d)*dt + Zd(:,d);
d = d+1;
end
% must renormalize
buf = Zd(:,1).'*TCdyn_single2(0,myorb(:,1));
Zd = Zd*1/buf;
Zd(:,1)./Zd(:,length(Zd))
PRC = Zd;
clear Zd;





for jj = 1:length(fa)
%%Get IRC
X1 = [0 0 0 0].';
X1(fa(jj)) = 1;
YY =  V.'
IRCinit = (YY\X1);
Zd = IRCinit;
myorb = fliplr(perorb.');

d = 1;
for k = 1:length(myorb)-1
[k mm]
Zd(:,d+1) = numericaljacobian(@TCdyn_single2,myorb(:,k)).'*Zd(:,d)*dt - kappa(jj)*eye(length(init))*Zd(:,d)*dt + Zd(:,d);
d = d+1;
end
Zd(:,1)./Zd(:,length(Zd))
IRC{jj} = Zd;
clear Zd;
end


for k = 1:length(perorb)
    [H1(:,:,k),H2(:,:,k),H3(:,:,k),H4(:,:,k)] = symsderivs(myorb(:,k));
    k
end


kappas = kappa;
clear kappa;
lambdas = lambda;
clear lambda;




%%get B functions
%%%%%%%%%%%%%%%%%%%%%%%
for ppp = 1:length(fa)
%%get all ircs
clear Zd_u;
%%%%Evluate IRC2
init = [0 0 0 0].';
myorb = fliplr(perorb.');
ydiffuse = flipud(ydiff{ppp}).';
Zd2 = PRC;
kappa = 0;
lambda = lambdas(ppp);

for mm = 1:3
   
d = 1;
Zd_u = init;
J = zeros(length(init),length(init));
   
for k = 1:length(perorb)-1
Zd_u(:,d+1) =(Zd2(1,k)*H1(:,:,k) + Zd2(2,k)*H2(:,:,k) + Zd2(3,k)*H3(:,:,k) + Zd2(4,k)*H4(:,:,k) )*ydiffuse(:,k)*dt + numericaljacobian(@TCdyn_single2,myorb(:,k)).'*Zd_u(:,d)*dt -  kappa*eye(length(init))*Zd_u(:,d)*dt + Zd_u(:,d);
d = d+1;
end

eps = 1e-6;
for p = 1:length(init)
    pert = zeros(size(init));
    pert(p) = eps;
    Zd_p = init + pert;
    d = 1;
    for k = 1:length(perorb)-1
    
Zd_p(:,d+1) =(Zd2(1,k)*H1(:,:,k) + Zd2(2,k)*H2(:,:,k) + Zd2(3,k)*H3(:,:,k) + Zd2(4,k)*H4(:,:,k) )*ydiffuse(:,k)*dt  +   numericaljacobian(@TCdyn_single2,myorb(:,k)).'*Zd_p(:,d)*dt - kappa*eye(length(init))*Zd_p(:,d)*dt + Zd_p(:,d);
d = d+1;

    end
   
    J(:,p) = ( Zd_p(:,length(Zd_p)) - Zd_u(:,length(Zd_u)) ).'/eps;
end
   
L1 = 1/lambda;
[V1,D1] = eig(J);
qq = find(abs(diag(D1)-L1) < .02);

J = J-L1*eye(length(init));
mydiff = Zd_u(:,length(Zd_u))-L1*Zd_u(:,1);
dx = J\-mydiff; 

init = init + dx;
[dx init]
end

d = 1;
Zd_u = init;
for k = 1:length(perorb)-1
Zd_u(:,d+1) = (Zd2(1,k)*H1(:,:,k) + Zd2(2,k)*H2(:,:,k) + Zd2(3,k)*H3(:,:,k) + Zd2(4,k)*H4(:,:,k)  )*ydiffuse(:,k)*dt   +   numericaljacobian(@TCdyn_single2,myorb(:,k)).'*Zd_u(:,d)*dt - kappa*eye(length(init))*Zd_u(:,d)*dt + Zd_u(:,d);
d = d+1; 
end

buf = (exp(-(log(lambda)/T)*timeorb));
buf = ones(length(init),1)*buf.';
Zd_u = buf.*Zd_u;
B{ppp} = Zd_u;

end


%%%%%get C functions%%%%%%
for ppp = 1:length(fa)
    for qqq = 1:length(fa)

%%get all ircs
clear Zd_u;
%%%%Evluate IRC2
init = [0 0 0 0].';
myorb = fliplr(perorb.');
ydiffuse = flipud(ydiff{qqq}).';
Zd2 = IRC{ppp};
kappa = kappas(ppp);
lambda = lambdas(qqq)

for mm = 1:3
   
d = 1;
Zd_u = init;
J = zeros(length(init),length(init));
   
for k = 1:length(perorb)-1
Zd_u(:,d+1) =(Zd2(1,k)*H1(:,:,k) + Zd2(2,k)*H2(:,:,k) + Zd2(3,k)*H3(:,:,k) + Zd2(4,k)*H4(:,:,k)  )*ydiffuse(:,k)*dt + numericaljacobian(@TCdyn_single2,myorb(:,k)).'*Zd_u(:,d)*dt -  kappa*eye(length(init))*Zd_u(:,d)*dt + Zd_u(:,d);
d = d+1;
end

eps = 1e-6;
for p = 1:length(init)
    pert = zeros(size(init));
    pert(p) = eps;
    Zd_p = init + pert;
    d = 1;
    for k = 1:length(perorb)-1
    
Zd_p(:,d+1) =(Zd2(1,k)*H1(:,:,k) + Zd2(2,k)*H2(:,:,k) + Zd2(3,k)*H3(:,:,k) + Zd2(4,k)*H4(:,:,k) )*ydiffuse(:,k)*dt  +   numericaljacobian(@TCdyn_single2,myorb(:,k)).'*Zd_p(:,d)*dt - kappa*eye(length(init))*Zd_p(:,d)*dt + Zd_p(:,d);
d = d+1;

    end
   
    J(:,p) = ( Zd_p(:,length(Zd_p)) - Zd_u(:,length(Zd_u)) ).'/eps;
end
   
L1 = 1/lambda;

[V1,D1] = eig(J);
qq = find(abs(diag(D1)-L1) < .02);

J = J-L1*eye(length(init));
mydiff = Zd_u(:,length(Zd_u))-L1*Zd_u(:,1);


if ppp==qqq
J2 = [J;1 0 0 0];
mydiff2 = [mydiff;0];
dx = J2\-mydiff2;
else
   dx = J\-mydiff; 
end  


init = init + dx;

[dx init]
end

if ppp==qqq
   J22 = numericaljacobian(@TCdyn_single2,myorb(:,1));
   I22 = Zd2(:,1);
   F22 = TCdyn_single2(0,myorb(:,1));
   difval = ydiffuse(:,1);
  
   exactval = I22'*(kappa*eye(length(init))-J22)*difval;
   isval = F22.'*Zd_u(:,1);
   canchgval = F22.'*V1(:,qq);
   Bval = (exactval-isval)/canchgval;
   init = init + Bval*V1(:,qq);
end

d = 1;
Zd_u = init;
for k = 1:length(perorb)-1
Zd_u(:,d+1) = (Zd2(1,k)*H1(:,:,k) + Zd2(2,k)*H2(:,:,k) + Zd2(3,k)*H3(:,:,k) + Zd2(4,k)*H4(:,:,k) )*ydiffuse(:,k)*dt   +   numericaljacobian(@TCdyn_single2,myorb(:,k)).'*Zd_u(:,d)*dt - kappa*eye(length(init))*Zd_u(:,d)*dt + Zd_u(:,d);
d = d+1; 
end

buf = (exp(-(log(lambda)/T)*timeorb));
buf = ones(length(init),1)*buf.';
Zd_u = buf.*Zd_u;
C{ppp,qqq} = Zd_u;


end
end



ydiff{1} = flipud(ydiff{1});
for k = 1:length(fa)
   ydiffall{k} = ydiff{k}.';
   lambda = lambdas(k);
   buf = (exp(-(log(lambda)/T)*timeorb));
    buf = ones(length(init),1)*buf.';
    ydiffall{k} = buf.*ydiffall{k};
    ydiffall{k} = fliplr(ydiffall{k});
end



return


%%now flip everything
PRC = fliplr(PRC);
for k = 1:length(fa)
    IRC{k} = fliplr(IRC{k});
    B{k} = fliplr(B{k})
end
for k = 1:length(fa)
    for j = 1:length(fa)
        C{k,j} = fliplr(C{k,j});
    end
end
timeorb = flipud(timeorb);



%save reduced model parameters
save('reducedparameters_noselfcoupling_dtp0001_Ib3p75_newnormalize','omega','T','timeorb','PRC','IRC','V','C','kappas','perorb','ydiffall','B')








