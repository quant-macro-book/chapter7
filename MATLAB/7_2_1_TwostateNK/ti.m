function [yvec0 pvec0 rvec0] = ti(m)

% disp('')
% disp('-+- Solve a two-state model with time iteration -+-');

%% STEP 1(a): �O���b�h����
Gs = [m.sH; m.sL];
Ps = [1-m.pH m.pH;
    1-m.pL m.pL]; 

%% ��͓I��
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

%% STEP 1(b): ����֐��̏����l�𓖂Đ���
Ns = 2;
% ��͓I���������l�Ƃ���(1��̌J��Ԃ��Ŏ���)
% yvec0 = [yH; yL];
% pvec0 = [piH; piL];
% rvec0 = [rH; 0];
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
        s0 = Gs(is);

        % �Â�����֐�������Ғl(ye, pie)���v�Z
        ye = Ps(is,:)*yvec0;
        pie = Ps(is,:)*pvec0;

        % ���Ғl�����^�Ƃ��čœK��
        r0 = max(m.rstar + m.phi*pie, 0);
        y0 = ye - (r0 - pie - s0);
        p0 = m.kap*y0 + m.bet*pie;

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