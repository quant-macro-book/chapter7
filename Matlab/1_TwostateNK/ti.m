function [yvec0 pvec0 rvec0] = ti(m)

% disp('')
% disp('-+- Solve a two-state model with time iteration -+-');

%% STEP 1(a): グリッド生成
Gs = [m.sH; m.sL];
Ps = [1-m.pH m.pH;
    1-m.pL m.pL]; 

%% 解析的解
A = [-1+(1-m.pH) m.pH -(m.phi-1)*(1-m.pH) -(m.phi-1)*m.pH;
m.kap 0 -1+m.bet*(1-m.pH) m.bet*m.pH;
(1-m.pL) -1+m.pL (1-m.pL) m.pL;
0 m.kap m.bet*(1-m.pL) -1+m.bet*m.pL];
b = [m.rstar-m.sH;0;-m.sL;0];
x = A\b;
yH  = x(1);
yL  = x(2);
piH = x(3);
piL = x(4);
rH = m.rstar + m.phi*((1-m.pH)*piH + m.pH*piL);

%% STEP 1(b): 政策関数の初期値を当て推量
Ns = 2;
% 解析的解を初期値とする(1回の繰り返しで収束)
% yvec0 = [yH; yL];
% pvec0 = [piH; piL];
% rvec0 = [rH; 0];
% 適当な初期値
yvec0 = zeros(Ns,1);
pvec0 = zeros(Ns,1);
rvec0 = zeros(Ns,1);
yvec1 = zeros(Ns,1);
pvec1 = zeros(Ns,1);
rvec1 = zeros(Ns,1);

%% STEP 4: 政策関数を繰り返し計算
diff = 1e+4; % 政策関数の繰り返し誤差
iter = 1; % ループ・カウンター

while(diff > m.tol)

    for is = 1:Ns

        % ショックの値
        s0 = Gs(is);

        % 古い政策関数から期待値(ye, pie)を計算
        ye = Ps(is,:)*yvec0;
        pie = Ps(is,:)*pvec0;

        % 期待値を所与として最適化
        r0 = max(m.rstar + m.phi*pie, 0);
        y0 = ye - (r0 - pie - s0);
        p0 = m.kap*y0 + m.bet*pie;

        % 新しい政策関数を保存
        yvec1(is,1) = y0;
        pvec1(is,1) = p0;
        rvec1(is,1) = r0;

    end
    
    % 繰り返し計算誤差を確認
    ydiff = max(abs(yvec1-yvec0));
    pdiff = max(abs(pvec1-pvec0));
    rdiff = max(abs(rvec1-rvec0));
    diff = max([ydiff pdiff rdiff]);

    disp([iter diff]);

    % 政策関数をアップデート
    yvec0 = yvec1;
    pvec0 = pvec1;
    rvec0 = rvec1;
    
    iter = iter + 1;
    
end