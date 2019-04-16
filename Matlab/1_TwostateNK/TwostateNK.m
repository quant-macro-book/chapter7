clear all;

m.rstar = 0.75; % pH=0のときの、定常状態での名目金利の値
m.bet = 1/(1+m.rstar/100); % 割引率(オイラー方程式の定常状態より)
m.phi = 5.0;  % テイラー係数(注: 小さいとiL=0にならない)
m.pL = 0.75;  % 危機の継続確率
m.sH = m.rstar; % 状態Hでの自然利子率の値

% カリブレーション
% yLとpiLのターゲットにpH=0のときのモデルの値を合わせるように、sLとkapの値をセット

m.pH = 0.0; % 危機が起こる確率
x0 = [-2.0, 0.01]; % sLとkapの初期値

% yLとpiLのターゲット
yLtar = -7.0;
piLtar = -1.0/4;

% 最小化関数(Matlabの場合fminsearch)を用いる
x = fminsearch(@dist,x0,[],m.sH,m.pH,m.pL,m.bet,m.phi,m.rstar,yLtar,piLtar);

% カリブレートしたパラメータをセット
m.sL = x(1); % 状態Lでの自然利子率の値
m.kap = x(2); % フィリップス曲線の傾き

m.maxiter = 2000; % 繰り返し回数の最大値
m.tol = 1e-5; % 許容誤差

%%
tic;
[yvec0 pvec0 rvec0] = ti(m);
toc;

m.pH = 0.025
tic;
[yvec1 pvec1 rvec1] = ti(m);
toc;

%%
xvec = [0 1];

figure;
subplot(231);
plot(xvec,rvec0*4,'k*-','LineWidth',3.0);
hold on;
plot(xvec,rvec1*4,'k*--','LineWidth',3.0);
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
% saveas(gcf,'simplepf.eps','epsc2');