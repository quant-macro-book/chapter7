function f = dist(x0,sH,pH,pL,bet,phi,rstar,yLtar,piLtar)

sL = x0(1);
kap = x0(2);

% 解析的解
A = [-1+(1-pH) pH -(phi-1)*(1-pH) -(phi-1)*pH;
    kap 0 -1+bet*(1-pH) bet*pH;
    (1-pL) -1+pL (1-pL) pL;
    0 kap bet*(1-pL) -1+bet*pL];
b = [rstar-sH;0;-sL;0];
x = A\b;
%yH  = x(1);
yL  = x(2);
%piH = x(3);
piL = x(4);

% データの値とモデルの値の間の距離を最小化
f = (yL-yLtar)^2 + (piL-piLtar)^2;