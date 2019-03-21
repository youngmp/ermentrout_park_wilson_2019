function[H1,H2,H3,H4,H5] = symsderivs(input)

%  syms vth hth rth sth
% 

vth = input(1);
hth = input(2);
rth = input(3);
sth = input(4);


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
Ib = 4;
%thalamic cell currents
    IL=gL*(vth-EL);
    INa=gNa*(m_inf.^3).*hth.*(vth-ENa);
    IK=gK*((0.75*(1-hth)).^4).*(vth-EK);
    IT=gT*(p_inf.^2).*rth.*(vth-ET);
    

    
    %synaptic
    alpha = 3;
Vt = -20;
sigmat = .8;
beta = .2;
gsyn = 0.015*15;  
Esyn = -61;
%%%

 
Isyn = 0*gsyn*(sth)*(vth-Esyn);

%Differential Equations for cells thalamic
%     vth_dot= 1/Cm*(-IL-INa-IK-IT+Ib - Isyn); 
%     hth_dot=(hth_inf-hth)./tauh;
%     rth_dot=(rth_inf-rth)./taur;
%     sth_dot = alpha*(1-sth)*(1/(1+exp(-(vth-Vt)/sigmat)))-beta*sth;




% S1 = 1/Cm*(-IL-INa-IK-IT+Ib - Isyn); 
% H1 = [diff(diff(S1,vth),vth) diff(diff(S1,vth),hth) diff(diff(S1,vth),rth) diff(diff(S1,vth),sth)
%     diff(diff(S1,hth),vth) diff(diff(S1,hth),hth) diff(diff(S1,hth),rth) diff(diff(S1,hth),sth)
%     diff(diff(S1,rth),vth) diff(diff(S1,rth),hth) diff(diff(S1,rth),rth) diff(diff(S1,rth),sth)
%     diff(diff(S1,sth),vth) diff(diff(S1,sth),hth) diff(diff(S1,sth),rth) diff(diff(S1,sth),sth) ];


H1 = [ (500*rth*vth*exp(- (5*vth)/31 - 300/31))/(961*(exp(- (5*vth)/31 - 300/31) + 1)^3) - (200*rth*exp(- (5*vth)/31 - 300/31))/(31*(exp(- (5*vth)/31 - 300/31) + 1)^3) - (18*hth*exp(- vth/7 - 37/7))/(7*(exp(- vth/7 - 37/7) + 1)^4) - (1500*rth*vth*exp(- (10*vth)/31 - 600/31))/(961*(exp(- (5*vth)/31 - 300/31) + 1)^4) + (9*hth*exp(- vth/7 - 37/7)*(vth - 50))/(49*(exp(- vth/7 - 37/7) + 1)^4) - (36*hth*exp(- (2*vth)/7 - 74/7)*(vth - 50))/(49*(exp(- vth/7 - 37/7) + 1)^5), - 3/(exp(- vth/7 - 37/7) + 1)^3 - 15*((3*hth)/4 - 3/4)^3 - (9*exp(- vth/7 - 37/7)*(vth - 50))/(7*(exp(- vth/7 - 37/7) + 1)^4), - 10/(exp(- (5*vth)/31 - 300/31) + 1)^2 - (100*vth*exp(- (5*vth)/31 - 300/31))/(31*(exp(- (5*vth)/31 - 300/31) + 1)^3), 0
                                                                                                                                                                                                                                                                                                                                                  - 3/(exp(- vth/7 - 37/7) + 1)^3 - 15*((3*hth)/4 - 3/4)^3 - (9*exp(- vth/7 - 37/7)*(vth - 50))/(7*(exp(- vth/7 - 37/7) + 1)^4),                                                                                       -(135*((3*hth)/4 - 3/4)^2*(vth + 90))/4,                                                                                                                      0, 0
                                                                                                                                                                                                                                                                                                                                                         - 10/(exp(- (5*vth)/31 - 300/31) + 1)^2 - (100*vth*exp(- (5*vth)/31 - 300/31))/(31*(exp(- (5*vth)/31 - 300/31) + 1)^3),                                                                                                                             0,                                                                                                                      0, 0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              0,                                                                                                                             0,                                                                                                                      0, 0];
 

% S1 = (hth_inf-hth)./tauh;
% H1 = [diff(diff(S1,vth),vth) diff(diff(S1,vth),hth) diff(diff(S1,vth),rth) diff(diff(S1,vth),sth)
%     diff(diff(S1,hth),vth) diff(diff(S1,hth),hth) diff(diff(S1,hth),rth) diff(diff(S1,hth),sth)
%     diff(diff(S1,rth),vth) diff(diff(S1,rth),hth) diff(diff(S1,rth),rth) diff(diff(S1,rth),sth)
%     diff(diff(S1,sth),vth) diff(diff(S1,sth),hth) diff(diff(S1,sth),rth) diff(diff(S1,sth),sth) ];


H2 = [ (exp(vth/4 + 41/4)*((8*exp(- vth/18 - 23/9))/1125 - (4*exp(- vth/5 - 23/5))/(5*(exp(- vth/5 - 23/5) + 1)^2)))/(2*(exp(vth/4 + 41/4) + 1)^2) - (hth - 1/(exp(vth/4 + 41/4) + 1))*((4*exp(- vth/18 - 23/9))/10125 - (4*exp(- vth/5 - 23/5))/(25*(exp(- vth/5 - 23/5) + 1)^2) + (8*exp(- (2*vth)/5 - 46/5))/(25*(exp(- vth/5 - 23/5) + 1)^3)) + (exp(vth/2 + 41/2)*((16*exp(- vth/18 - 23/9))/125 + 4/(exp(- vth/5 - 23/5) + 1)))/(8*(exp(vth/4 + 41/4) + 1)^3) - (exp(vth/4 + 41/4)*((16*exp(- vth/18 - 23/9))/125 + 4/(exp(- vth/5 - 23/5) + 1)))/(16*(exp(vth/4 + 41/4) + 1)^2), (8*exp(- vth/18 - 23/9))/1125 - (4*exp(- vth/5 - 23/5))/(5*(exp(- vth/5 - 23/5) + 1)^2), 0, 0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         (8*exp(- vth/18 - 23/9))/1125 - (4*exp(- vth/5 - 23/5))/(5*(exp(- vth/5 - 23/5) + 1)^2),                                                                                       0, 0, 0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0,                                                                                       0, 0, 0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0,                                                                                       0, 0, 0];
 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
% S1 = (rth_inf-rth)./taur;
% H1 = [diff(diff(S1,vth),vth) diff(diff(S1,vth),hth) diff(diff(S1,vth),rth) diff(diff(S1,vth),sth)
%     diff(diff(S1,hth),vth) diff(diff(S1,hth),hth) diff(diff(S1,hth),rth) diff(diff(S1,hth),sth)
%     diff(diff(S1,rth),vth) diff(diff(S1,rth),hth) diff(diff(S1,rth),rth) diff(diff(S1,rth),sth)
%     diff(diff(S1,sth),vth) diff(diff(S1,sth),hth) diff(diff(S1,sth),rth) diff(diff(S1,sth),sth) ];


H3 =  [ exp(vth/2 + 42)/(8*(exp(vth/4 + 21) + 1)^3*(exp(- vth/20 - 5/4) + 28)) - exp(vth/4 + 21)/(16*(exp(vth/4 + 21) + 1)^2*(exp(- vth/20 - 5/4) + 28)) - (exp(- vth/10 - 5/2)*(rth - 1/(exp(vth/4 + 21) + 1)))/(200*(exp(- vth/20 - 5/4) + 28)^3) + (exp(- vth/20 - 5/4)*(rth - 1/(exp(vth/4 + 21) + 1)))/(400*(exp(- vth/20 - 5/4) + 28)^2) - (exp(vth/4 + 21)*exp(- vth/20 - 5/4))/(40*(exp(vth/4 + 21) + 1)^2*(exp(- vth/20 - 5/4) + 28)^2), 0, -exp(- vth/20 - 5/4)/(20*(exp(- vth/20 - 5/4) + 28)^2), 0
                                                                                                                                                                                                                                                                                                                                                                                                                                        0, 0,                                                      0, 0
                                                                                                                                                                                                                                                                                                                                                                                   -exp(- vth/20 - 5/4)/(20*(exp(- vth/20 - 5/4) + 28)^2), 0,                                                      0, 0
                                                                                                                                                                                                                                                                                                                                                                                                                                        0, 0,                                                      0, 0];
 


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
% S1 = alpha*(1-sth)*(1/(1+exp(-(vth-Vt)/sigmat)))-beta*sth;
% H1 = [diff(diff(S1,vth),vth) diff(diff(S1,vth),hth) diff(diff(S1,vth),rth) diff(diff(S1,vth),sth)
%     diff(diff(S1,hth),vth) diff(diff(S1,hth),hth) diff(diff(S1,hth),rth) diff(diff(S1,hth),sth)
%     diff(diff(S1,rth),vth) diff(diff(S1,rth),hth) diff(diff(S1,rth),rth) diff(diff(S1,rth),sth)
%     diff(diff(S1,sth),vth) diff(diff(S1,sth),hth) diff(diff(S1,sth),rth) diff(diff(S1,sth),sth) ];
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
H4 =     [ (25*exp(- (5*vth)/4 - 25)*(3*sth - 3))/(16*(exp(- (5*vth)/4 - 25) + 1)^2) - (25*exp(- (5*vth)/2 - 50)*(3*sth - 3))/(8*(exp(- (5*vth)/4 - 25) + 1)^3), 0, 0, -(15*exp(- (5*vth)/4 - 25))/(4*(exp(- (5*vth)/4 - 25) + 1)^2)
                                                                                                                                                    0, 0, 0,                                                             0
                                                                                                                                                    0, 0, 0,                                                             0
                                                                                        -(15*exp(- (5*vth)/4 - 25))/(4*(exp(- (5*vth)/4 - 25) + 1)^2), 0, 0,                                                             0];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
