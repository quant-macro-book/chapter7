clear all;

bet = 0.99;
kap = 0.03;
rstar = 100/bet-100;
% no labor: govt mp is always equal to one!!!
% N = 0.0;
% gam = 0.0;
N = 1/3;
gam = 0.29;
phi = 1.5;
gy = 0.2;
sig = 2;
rho = 0.80;

zlbflag = 1;

nrho = 121;
% rhovec = linspace(0.0,0.11,nrho)';
rhovec = linspace(0.7,0.81,nrho)';
mpvec0 = zeros(nrho,2);
mpvec1 = zeros(nrho,2);

for i = 1:nrho
    rho = rhovec(i);
    [yL1 pL1] = cer_nl_v2(1.0,bet,kap,rstar,N,gam,phi,gy,sig,rho,0);
    [yL0 pL0] = cer_nl_v2(0.0,bet,kap,rstar,N,gam,phi,gy,sig,rho,0);
    mpvec0(i,1) = (yL1-yL0)/gy;
    mpvec0(i,2) = (pL1-pL0);
    [yL1 pL1] = cer_nl_v2(1.0,bet,kap,rstar,N,gam,phi,gy,sig,rho,1);
    [yL0 pL1] = cer_nl_v2(0.0,bet,kap,rstar,N,gam,phi,gy,sig,rho,1);
    mpvec1(i,1) = (yL1-yL0)/gy;
    mpvec1(i,2) = (pL1-pL0);
end

figure;
plot(rhovec,mpvec0(:,2),'r*--');
% hold on;
% plot(rhovec,mpvec1(:,2),'b*-');
% legend('without ZLB','with ZLB');
xlim([rhovec(1) rhovec(nrho)]);
% ylim([1.0 7.0]);
