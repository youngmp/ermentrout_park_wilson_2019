function dy = isostablereducedmodel(t,y)

global T1 T2 F1 F2 F3 F4 F5 F6 alpha kappas


theta1 = mod(y(1),2*pi);
psi1 = y(2);
theta2 = mod(y(3),2*pi);
psi2 = y(4);

dy = zeros(4,1);


dy(1) = psi1 * alpha * qinterp2(T1,T2,F2,theta2,theta1) + psi2 * alpha * qinterp2(T1,T2,F3,theta2,theta1);
dy(2) = kappas*psi1 + alpha * qinterp2(T1,T2,F4,theta2,theta1) + alpha * psi1 * qinterp2(T1,T2,F5,theta2,theta1) + alpha * psi2 * qinterp2(T1,T2,F6,theta2,theta1);

dy(3) = alpha * psi2 * qinterp2(T1,T2,F2,theta1,theta2) + alpha * psi1 * qinterp2(T1,T2,F3,theta1,theta2);
dy(4) = kappas*psi2 + alpha * qinterp2(T1,T2,F4,theta1,theta2) + alpha * psi2 * qinterp2(T1,T2,F5,theta1,theta2) + alpha * psi1 * qinterp2(T1,T2,F6,theta1,theta2);


