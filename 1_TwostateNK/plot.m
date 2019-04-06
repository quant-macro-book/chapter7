clear all;

%bet = 0.9950;
%kap = 0.024;
%phi = 2.0;
%pibar = 1.005;
%rstar = log(pibar/bet)*100;

load pfmatpH0.025.mat yvec0 pvec0 ivec0 Gs;
yvec1 = yvec0;
pvec1 = pvec0;
ivec1 = ivec0;
% yvec1 = yvec0([2 1]);
% pvec1 = pvec0([2 1]);
% ivec1 = ivec0([2 1]);
load pfmatpH0.0.mat yvec0 pvec0 ivec0 Gs;
% yvec0 = yvec0([2 1]);
% pvec0 = pvec0([2 1]);
% ivec0 = ivec0([2 1]);
xvec = [0 1];

figure;
subplot(231);
plot(xvec,ivec0*4,'k*-','LineWidth',3.0);
hold on;
plot(xvec,ivec1*4,'k*--','LineWidth',3.0);
plot(xvec,[0 0],'r-');
title('政策金利');
xticks([0 1]);
xticklabels({'H','L'});
set(gca,'Fontsize',12);

subplot(232);
plot(xvec,yvec0,'k*-','LineWidth',3.0);
hold on;
plot(xvec,yvec1,'k*--','LineWidth',3.0);
plot(xvec,[0 0],'r-');
title('産出ギャップ');
xticks([0 1]);
xticklabels({'H','L'});
set(gca,'Fontsize',12);

subplot(233);
plot(xvec,pvec0*4,'k*-','LineWidth',3.0);
hold on;
plot(xvec,pvec1*4,'k*--','LineWidth',3.0);
plot(xvec,[0 0],'r-');
title('インフレ率');
xticks([0 1]);
xticklabels({'H','L'});
set(gca,'Fontsize',12);
m = legend('p_H=0','p_H=0.025','Location','SouthWest');
m.FontSize = 8;
saveas (gcf,'simplepf.eps','epsc2');