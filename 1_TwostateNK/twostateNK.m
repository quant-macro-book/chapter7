%% ショックの状態が2つのときの準線形モデルにおける、時間反復法(time iteration method)の解法
clear all;

%% *** カリブレーション ***
% todo: juliaファイルのように構造体に書き換える？
rstar = 0.75; % 状態Hのときの自然利子率の値; pH=0のときは名目金利の値とも等しい
bet = 1/(1+rstar/100); % 割引率(オイラー方程式の定常状態より)
phi = 5.0; % テイラー係数(注: 小さいとiL=0にならない)
pL = 0.75; % 危機の継続確率
sH = rstar;

% カリブレーション
% yLとpiLのターゲットにpH=0のときのモデルの値を合わせるように、sLとkapの値をセット
% 最小化関数fminsearchを用いる
pH = 0.0; % 危機が起こる確率
x0 = [-2.0; 0.01]; % sLとkapの初期値
yLtar = -7.0; %[-7.0; -1.0/4]; % yLとpiLのターゲット
piLtar = -1.0/4;

[x fval flag] = fminsearch(@dist,x0,[],sH,pH,pL,bet,phi,rstar,yLtar,piLtar);
sL = x(1); % 状態Lのときの自然利子率の値
kap = x(2); % フィリップス曲線の傾き
% for check
sL
kap
dist(x,sH,pH,pL,bet,phi,rstar,yLtar,piLtar)
pause;

%pH = 0.025;
%=================

%% 解析的解
A = [-1+(1-pH) pH -(phi-1)*(1-pH) -(phi-1)*pH;
    kap 0 -1+bet*(1-pH) bet*pH;
    (1-pL) -1+pL (1-pL) pL;
    0 kap bet*(1-pL) -1+bet*pL];
b = [rstar-sH;0;-sL;0];
x = A\b;
yH  = x(1);
yL  = x(2);
piH = x(3);
piL = x(4);
% check>0
iH = rstar + phi*((1-pH)*piH + pH*piL);

%% 数値的解の計算開始

tic

disp('')
disp('-+- Solve a two-state model with time iteration -+-');

%% STEP 1(a): グリッド生成
Gs = [sH; sL];
Ps = [1-pH pH; 1-pL pL]; 

%% STEP 1(b): 政策関数の初期値を当て推量
Ns = 2;
% 解析的解を初期値とする(1回の繰り返しで収束)
yvec0 = [yH; yL];
pvec0 = [piH; piL];
ivec0 = [iH; 0];
% 適当な初期値
% yvec0 = zeros(Ns,1);
% pvec0 = zeros(Ns,1);
% ivec0 = zeros(Ns,1);
yvec1 = zeros(Ns,1);
pvec1 = zeros(Ns,1);
ivec1 = zeros(Ns,1);

% *** 収束の基準 ***
crit = 1e-5; % 許容誤差(STEP 2)
diff = 1e+4; % 政策関数の繰り返し誤差
iter = 1 % ループ・カウンター

%% STEP 4: 政策関数を繰り返し計算

while(diff>crit)

    for is = 1:Ns

        % ショックの値
        s0 = Gs(is);

        % 古い政策関数から期待値(ye, pie)を計算
        ye = Ps(is,:)*yvec0;
        pie = Ps(is,:)*pvec0;

        % 期待値を所与として最適化
        i0 = max(0, rstar + phi*pie);
        y0 = ye - (i0 - pie - s0);
        p0 = kap*y0 + bet*pie;

        % 新しい政策関数を保存
        yvec1(is,1) = y0;
        pvec1(is,1) = p0;
        ivec1(is,1) = i0;

    end
    
    % 繰り返し計算誤差を確認
    ydiff = max(abs(yvec1-yvec0));
    pdiff = max(abs(pvec1-pvec0));
    idiff = max(abs(ivec1-ivec0));
    diff = max([ydiff pdiff idiff]);

    disp([iter diff]);

    % 政策関数をアップデート
    yvec0 = yvec1;
    pvec0 = pvec1;
    ivec0 = ivec1;
    
    iter = iter + 1;
    
end

% %% 計算結果をコマンドウィンドウに表示
% 
% disp('-+- PARAMETER VALUES -+-');
% disp('');
% fprintf('beta=%5.2f, kappa=%5.2f, phi=%5.2f, rstar=%5.2f \n', bet, kap, phi, rstar);
disp('');