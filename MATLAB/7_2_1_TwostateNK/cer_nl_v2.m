function [yL pL] = cer_nl_v2(gL,bet,kap,rstar,N,gam,phi,gy,sig,rho,zlbflag)

fprintf('gL = %1.4f\n',gL);

dc = gam*(1-sig)-1;
dn = -(1-gam)*(1-sig)*N/(1-N);

ng = 2;
pH = 0.0;
pL = rho;
sH = rstar;
sL = -5*1.5; % irrelevant for govt mp as the ZLB is binding in the low state
Gg = [0.0; gL];
Pg = [1-pH pH; 1-pL pL];
Gs = [sH; sL];

lamvec0 = zeros(ng,1);
pivec0  = zeros(ng,1);
lamvec1 = zeros(ng,1);
pivec1  = zeros(ng,1);
cvec = zeros(ng,1);
nvec = zeros(ng,1);
rvec = zeros(ng,1);

iter = 0;
diff = 10000;

while (diff>1e-10)
    
    for ig = 1:ng

        g0 = Gg(ig);
        lame = Pg(ig,:)*lamvec0;
        pie = Pg(ig,:)*pivec0;

        if (zlbflag)
            r0 = max(rstar + phi*pie, 0);
        else
            r0 = rstar + phi*pie;
        end
        lam0 = lame + r0 - pie - Gs(ig); 

        c0 = (lam0 - dn*gy*g0)/(dc + dn*(1-gy));
        n0 = (1-gy)*c0 + gy*g0;
        pi0 = bet*pie + kap*(c0 + N/(1-N)*n0);

        lamvec1(ig) = lam0;
        pivec1(ig) = pi0;
        cvec(ig) = c0;
        nvec(ig) = n0;
        rvec(ig) = r0;

    end
    
    diff_lam = max(abs(lamvec1-lamvec0));
    diff_pi  = max(abs(pivec1-pivec0));
    diff = max(diff_lam,diff_pi);    
    iter = iter + 1;

    if (mod(iter,1000)==0) disp([iter diff_lam diff_pi]); end
    
    lamvec0 = lamvec1;
    pivec0 = pivec1;
    
end

if (rvec(2)>0); disp('ZLB is not binding at state L!!!'); end;

yL = nvec(2);
pL = pivec1(2);