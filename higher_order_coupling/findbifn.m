clear all
close all
clc

global gsyn

init = zeros(4*2,1) + [.01 0 0 0 0 0 0 0].';
options = odeset('abstol',1e-8,'reltol',1e-8);
tspan = [0 20000];

maxg = 0.0463;
ming = 0.0462;

gsyn = (maxg+ming)/2;
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

plot(t,Y(:,1))
hold on
plot(t,Y(:,5))
getframe

gsyn

return
