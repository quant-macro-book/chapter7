% Replication of Adam and Billi (2007,JME)
% January 2018, Takeki Sunakawa
clear all;

m.rstar = 3.5/4; % pH=0�̂Ƃ��́A����Ԃł̖��ڋ����̒l
m.bet = 1/(1+m.rstar/100); % ������(�I�C���[�������̒���Ԃ��)
m.sig = 6.0; %6.25;
m.alp = 0.66;
m.the = 7.66;
m.ome = 0.47;
m.kap = (1-m.alp)*(1-m.alp*m.bet)/m.alp*(1/m.sig+m.ome)/(1+m.ome*m.the); %0.024;
m.lam = 0.048/16; %0.003;

% joint shock process
m.rhou = 0;
m.rhog = 0.8;
m.sigu = 1e-5; %0.154;
m.sigg = 1.524;
m.Ng = 31;
m.Nu = 31;

m.maxiter = 2000; % �J��Ԃ��񐔂̍ő�l
m.tol = 1e-5; % ���e�덷

%%
tic;
[yvec0 pvec0 rvec0 Gg Gu] = ti(m);
toc;

%%
% replicate Figures 4-5 in the paper
for ig = 1:m.Ng
    
    for iu = 1:m.Nu
        
        ymat0(ig,iu) = yvec0(m.Nu*(ig-1)+iu,1);
        pmat0(ig,iu) = pvec0(m.Nu*(ig-1)+iu,1);
        rmat0(ig,iu) = rvec0(m.Nu*(ig-1)+iu,1);
        
    end
    
end

rmat0 = rmat0; %-rstar;
Gg = 4*(1/m.sig)*Gg;
%Gg = Gg-m.sig*m.rstar;
idu = ceil(m.Nu/2);
stg = 1;
edg = ceil(m.Ng/2);

figure;
subplot(311);
plot(Gg(stg:edg),ymat0(stg:edg,idu),'k-','LineWidth',2.0);
xlim([Gg(stg) Gg(edg)]); xticks([-1:1:3]);
ylim([-8 2]); yticks([-8:2:2]);
grid on;
ylabel('産出ギャップ, y');
subplot(312);
plot(Gg(stg:edg),4*pmat0(stg:edg,idu),'k-','LineWidth',2.0);
xlim([Gg(stg) Gg(edg)]); xticks([-1:1:3]);
ylim([-2 0.5]); yticks([-2:0.5:0.5]);
grid on;
ylabel('インフレ率, \pi');
subplot(313);
plot(Gg(stg:edg),4*rmat0(stg:edg,idu),'k-','LineWidth',2.0);
hold on;
plot(Gg(stg:edg),Gg(stg:edg),'k:');
xlim([Gg(stg) Gg(edg)]); xticks([-1:1:3]);
%ylim([-4 0]); yticks([-4:1:0]);
grid on;
xlabel('自然利子率, g');
ylabel('政策金利, r_{n}');