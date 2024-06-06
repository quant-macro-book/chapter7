function [yvec0 pvec0 rvec0 Gg Gu] = ti(m)

% disp('')
% disp('-+- Solve a two-state model with time iteration -+-');

%% STEP 1(a): グリッド生成
[Gg Pg] = tauchen(m.Ng,m.sig*m.rstar,m.rhog,m.sigg,3.0);
[Gu Pu] = tauchen(m.Nu,0,m.rhou,m.sigu,3.0);

Ns = m.Ng*m.Nu;
Gs = zeros(Ns,2);

for ig = 1:m.Ng
    
    for iu = 1:m.Nu
        
        Gs(m.Nu*(ig-1)+iu,1) = Gg(ig);
        Gs(m.Nu*(ig-1)+iu,2) = Gu(iu);
        
    end
    
end

Ps = kron(Pg,Pu);

%% STEP 1(b): 政策関数の初期値を当て推量
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
        g0 = Gs(is,1);
        u0 = Gs(is,2);

        % 古い政策関数から期待値(ye, pie)を計算
        ye = Ps(is,:)*yvec0;
        pie = Ps(is,:)*pvec0;

        % 期待値を所与として最適化
        p0 = (m.bet*pie+u0)/(1+m.kap^2/m.lam);
        y0 = (-m.kap/m.lam)*p0;
        r0 = (1/m.sig)*(ye - y0 + g0) + pie;

        if (r0<0)
            y0 = ye - m.sig*(0 - pie) + g0;
            p0 = m.kap*y0 + m.bet*pie + u0;
            r0 = 0;
        end

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