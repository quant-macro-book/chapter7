function [yvec0 pvec0 rvec0 Gg Gu] = ti(m)

% disp('')
% disp('-+- Solve a two-state model with time iteration -+-');

%% STEP 1(a): �O���b�h����
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

%% STEP 1(b): ����֐��̏����l�𓖂Đ���
% �K���ȏ����l
yvec0 = zeros(Ns,1);
pvec0 = zeros(Ns,1);
rvec0 = zeros(Ns,1);
yvec1 = zeros(Ns,1);
pvec1 = zeros(Ns,1);
rvec1 = zeros(Ns,1);

%% STEP 4: ����֐����J��Ԃ��v�Z
diff = 1e+4; % ����֐��̌J��Ԃ��덷
iter = 1; % ���[�v�E�J�E���^�[

while(diff > m.tol)

    for is = 1:Ns

        % �V���b�N�̒l
        g0 = Gs(is,1);
        u0 = Gs(is,2);

        % �Â�����֐�������Ғl(ye, pie)���v�Z
        ye = Ps(is,:)*yvec0;
        pie = Ps(is,:)*pvec0;

        % ���Ғl�����^�Ƃ��čœK��
        p0 = (m.bet*pie+u0)/(1+m.kap^2/m.lam);
        y0 = (-m.kap/m.lam)*p0;
        r0 = (1/m.sig)*(ye - y0 + g0) + pie;

        if (r0<0)
            y0 = ye - m.sig*(0 - pie) + g0;
            p0 = m.kap*y0 + m.bet*pie + u0;
            r0 = 0;
        end

        % �V��������֐���ۑ�
        yvec1(is,1) = y0;
        pvec1(is,1) = p0;
        rvec1(is,1) = r0;

    end
    
    % �J��Ԃ��v�Z�덷���m�F
    ydiff = max(abs(yvec1-yvec0));
    pdiff = max(abs(pvec1-pvec0));
    rdiff = max(abs(rvec1-rvec0));
    diff = max([ydiff pdiff rdiff]);

    disp([iter diff]);

    % ����֐����A�b�v�f�[�g
    yvec0 = yvec1;
    pvec0 = pvec1;
    rvec0 = rvec1;
    
    iter = iter + 1;
    
end