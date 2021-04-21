function x = Sigmaafun(nus, nua, Sigmasur,beta,alpha,delta,k,t,r)
%%% determines stress in the active material
vsp = 1+nus;
vsn = 1-nus;
vs2 = 1-2*nus;

vap = 1+nua;
van = 1-nua;
va2 = 1-2*nua;

D1 = 1+delta;

Cint1 = Cintfun(1,t,k);
CintD = Cintfun(D1,t,k);

Aa = va2 * (((-(van * alpha * (vs2 + 1) * beta) / 0.3e1 + (vsn * (D1 ^ 2 + vs2)) / 0.3e1) * vsp...
    - (beta * vap * vsn * (D1 - 1) * (D1 + 1)) / 0.3e1) * Cint1 + ...
    van * (vs2 + 1) * vsp * ((D1 ^ 2 * Sigmasur * vsn) + (CintD * beta * alpha) / 0.3e1))...
    * vap / van / ((vap * va2 * (D1 - 1) * (D1 + 1) * beta) + vsp * (D1 ^ 2 + vs2)) / vsn;

x = Aa/va2/vap-Cintfun(r,t,k)/3/r^2/van;

end

