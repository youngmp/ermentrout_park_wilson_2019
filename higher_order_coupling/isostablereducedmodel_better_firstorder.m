function dy = isostablereducedmodel_better(t,y)

global  alpha kappas F1_t1mt2 F2_t1mt2 F3_t1mt2 F4_t1mt2 F5_t1mt2 F6_t1mt2

TH = linspace(-2*pi,2*pi,length(F1_t1mt2));
PHI = y(1);
psi1 = y(2);
psi2 = y(3);

dy = zeros(3,1);

% dy(1) = alpha*interp1q(TH(:),F1_t1mt2(:),PHI) -  alpha*interp1q(TH(:),F1_t1mt2(:),-PHI);
% 
% dy(2) = kappas*psi1 +  alpha*interp1q(TH(:),F4_t1mt2(:),PHI);
% dy(3) = kappas*psi2 +  alpha*interp1q(TH(:),F4_t1mt2(:),-PHI);



dy(1) = alpha*interp1q(TH(:),F1_t1mt2(:),PHI) -  alpha*interp1q(TH(:),F1_t1mt2(:),-PHI);
dy(2) = kappas*psi1 +  alpha*interp1q(TH(:),F4_t1mt2(:),PHI) ;
dy(3) = kappas*psi2 +  alpha*interp1q(TH(:),F4_t1mt2(:),-PHI);


% dy(1) = psi1 * alpha * qinterp2(T1,T2,F2,theta2,theta1) + psi2 * alpha * qinterp2(T1,T2,F3,theta2,theta1);
% dy(2) = kappas*psi1 + alpha * qinterp2(T1,T2,F4,theta2,theta1) + alpha * psi1 * qinterp2(T1,T2,F5,theta2,theta1) + alpha * psi2 * qinterp2(T1,T2,F6,theta2,theta1);
% dy(3) = alpha * psi2 * qinterp2(T1,T2,F2,theta1,theta2) + alpha * psi1 * qinterp2(T1,T2,F3,theta1,theta2);
% dy(4) = kappas*psi2 + alpha * qinterp2(T1,T2,F4,theta1,theta2) + alpha * psi2 * qinterp2(T1,T2,F5,theta1,theta2) + alpha * psi1 * qinterp2(T1,T2,F6,theta1,theta2);


