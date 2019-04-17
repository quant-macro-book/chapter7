% replicate Figures 4-5 in the paper
for ig = 1:Ng
    
    for iu = 1:Nu
        
        ymat0(ig,iu) = yvec0(Nu*(ig-1)+iu,1);
        pmat0(ig,iu) = pvec0(Nu*(ig-1)+iu,1);
        imat0(ig,iu) = ivec0(Nu*(ig-1)+iu,1);
        
    end
    
end

imat0 = imat0; %-rstar;
Gg = Gg; %-sig*rstar;
idu = ceil(Nu/2);
stg = 1;
edg = ceil(Ng/2);

figure;
subplot(311);
plot(Gg(stg:edg),ymat0(stg:edg,idu));
xlim([Gg(stg) Gg(edg)]);
ylim([-8 2]); yticks([-8:2:2]);
grid on;
xlabel('g'); ylabel('y');
subplot(312);
plot(Gg(stg:edg),4*pmat0(stg:edg,idu));
xlim([Gg(stg) Gg(edg)]);
ylim([-2 0.5]); yticks([-2:0.5:0.5]);
grid on;
xlabel('g'); ylabel('\pi');
subplot(313);
plot(Gg(stg:edg),4*imat0(stg:edg,idu));
xlim([Gg(stg) Gg(edg)]);
%ylim([-4 0]); yticks([-4:1:0]);
grid on;
xlabel('g'); ylabel('i');

% figure;
% subplot(311);
% plot(Gu,ymat0(edg,:));
% xlim([Gu(1) Gu(end)]);
% ylim([-4 4]); yticks([-4:2:4]);
% grid on;
% xlabel('u'); ylabel('y');
% subplot(312);
% plot(Gu,4*pmat0(edg,:));
% xlim([Gu(1) Gu(end)]);
% ylim([-2 2]); yticks([-2:1:2]);
% grid on;
% xlabel('u'); ylabel('\pi');
% subplot(313);
% plot(Gu,4*imat0(edg,:));
% xlim([Gu(1) Gu(end)]);
% ylim([-3 2]); yticks([-3:1:2]);
% grid on;
% xlabel('u'); ylabel('i');
