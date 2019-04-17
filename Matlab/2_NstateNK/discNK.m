% Replication of Adam and Billi (2007,JME)
% January 2018, Takeki Sunakawa
clear all;

damp = 1.0;

% from Table 1 in the paper
rstar = 3.5/4;
bet = 1/(1+rstar/100);
sig = 6.25;
alp = 0.66;
the = 7.66;
ome = 0.47;
kap = (1-alp)*(1-alp*bet)/alp*(1/sig+ome)/(1+ome*the); %0.024;
lam = 0.048/16; %0.003;

% joint shock process
rhou = 0;
rhog = 0.8;
sigu = 0.154;
sigg = 1.524;
Ng = 31;
Nu = 31;
[Gg Pg] = tauchen(Ng,sig*rstar,rhog,sigg,3.0);
[Gu Pu] = tauchen(Nu,0,rhou,sigu,3.0);

Ns = Ng*Nu;
Gs = zeros(Ns,2);

for ig = 1:Ng
    
    for iu = 1:Nu
        
        Gs(Nu*(ig-1)+iu,1) = Gg(ig);
        Gs(Nu*(ig-1)+iu,2) = Gu(iu);
        
    end
    
end

Ps = kron(Pg,Pu);

% policy function iteration
yvec0 = zeros(Ns,1);
pvec0 = zeros(Ns,1);
rvec0 = zeros(Ns,1);
yvec1 = zeros(Ns,1);
pvec1 = zeros(Ns,1);
rvec1 = zeros(Ns,1);

tol = 1e-10;
diff = 1e+4;
iter = 0;

while(diff>tol)

    for is = 1:Ns

        g0 = Gs(is,1);
        u0 = Gs(is,2);

        ey = Ps(is,:)*yvec0;
        epi = Ps(is,:)*pvec0;

        p0 = (bet*epi+u0)/(1+kap^2/lam);
        y0 = (-kap/lam)*p0;
        r0 = (1/sig)*(ey - y0 + g0) + epi;
                
        if (r0<0)

            y0 = ey - sig*(0 - epi) + g0;
            p0 = kap*y0 + bet*epi + u0;
            r0 = 0;

        end
        
        yvec1(is,1) = y0;
        pvec1(is,1) = p0;
        rvec1(is,1) = r0;

    end
    
    ydiff = max(abs(yvec1-yvec0));
    pdiff = max(abs(pvec1-pvec0));
    rdiff = max(abs(rvec1-rvec0));
    diff = max([ydiff pdiff rdiff]);
    iter = iter + 1;
    s = sprintf( ' iteration %4d:  (%5.10f, %5.10f, %5.10f)', ...
        iter, ydiff, pdiff, rdiff);    
    disp(s);
    
    yvec0 = damp*yvec1 + (1-damp)*yvec0;
    pvec0 = damp*pvec1 + (1-damp)*pvec0;
    rvec0 = damp*rvec1 + (1-damp)*rvec0;
    
end

%save ab5.mat Gs yvec0 pvec0 ivec0;
%plotdiscab;