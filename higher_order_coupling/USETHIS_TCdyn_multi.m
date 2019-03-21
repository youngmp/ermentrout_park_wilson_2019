function [dy] = TCdyn_multi(t,Y)

global gsyn

dy = zeros(size(Y));

N = length(Y)/4;
for j = 1:N

vth = Y(1 + 4*(j-1))*100;
hth = Y(2 + 4*(j-1));
rth = Y(3 + 4*(j-1));
sth = Y(4 + 4*(j-1));

%thalamic cell parameters
    Cm = 1;
    gL = 0.15;
    EL = -75;
    gNa = 3;
    ENa = 50;
    gK = 5;
    EK = -90;
    gT = 10;
    ET = 0;
    hth_inf = 1/(1+exp((vth+41)/4));
    rth_inf = 1/(1+exp((vth+84)/4));
    ah = 0.128*exp(-(vth+46)/18);
    bh = 4/(1+exp(-(vth+23)/5));
    tauh = 1/(ah+bh);
    taur = 1*(28+exp(-(vth+25)/20));
    m_inf = 1/(1+exp(-(vth+37)/7));
    p_inf = 1/(1+exp(-(vth+60)/6.2));
Ib = 3.75;
%thalamic cell currents
    IL=gL*(vth-EL);
    INa=gNa*(m_inf.^3).*hth.*(vth-ENa);
    IK=gK*((0.75*(1-hth)).^4).*(vth-EK);
    IT=gT*(p_inf.^2).*rth.*(vth-ET);
    

%synaptic
    alpha = 3;
Vt = -20;
sigmat = .8;
beta = .20;
% gsyn = 0.1;   
Esyn = -60;
%%%
    
%     %synaptic
%     alpha = 3;
% Vt = -20;
% sigmat = .8;
% beta = .2;
% gsyn = 0.08;   
% Esyn = -61;
% %%%




    Isyn = (  gsyn*sum(Y(4:4:4*N))*(vth-Esyn)  -   0*gsyn*sth*(vth-Esyn)        )   ;

%Differential Equations for cells thalamic
    vth_dot= 1/Cm*(-IL-INa-IK-IT+Ib-Isyn)/100;   %- .01*(vth-sum(Y(1:4:4*N)));
    hth_dot=(hth_inf-hth)./tauh;
    rth_dot=(rth_inf-rth)./taur;
    sth_dot = alpha*(1-sth)*(1/(1+exp(-(vth-Vt)/sigmat)))-beta*sth;
    
    dy(1 + 4*(j-1)) = vth_dot;
    dy(2 + 4*(j-1)) = hth_dot;
    dy(3 + 4*(j-1)) = rth_dot;
    dy(4 + 4*(j-1)) = sth_dot;


end

dy = dy(:);
    