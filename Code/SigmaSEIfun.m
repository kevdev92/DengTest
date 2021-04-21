function x = SigmaSEIfun(nus, nua, Sigmasur,beta,alpha,delta,k,t,r)
%%% determines stress in the SEI

vsp = 1+nus;
vsn = 1-nus;
vs2 = 1-2*nus;

vap = 1+nua;
van = 1-nua;
va2 = 1-2*nua;

D1 = 1+delta;

Cint1 = Cintfun(1,t,k);
CintD = Cintfun(D1,t,k);

ASEI = vs2 * vsp * (-((alpha * beta * va2 * van * vap + vsp * alpha * van - vap *...
    vsn * (va2 + 1)) * beta * Cint1) / 0.3e1 + van * ((D1 ^ 2 * Sigmasur * vsn)...
    + (CintD * beta * alpha) / 0.3e1) * (beta * va2 * vap + vsp)) / van / ...
    (vap * va2 * (D1 - 1) * (D1 + 1) * beta + vsp * (D1 ^ 2 + vs2)) / vsn / beta;


BS =  vsp * (-((alpha * beta * va2 * van * vap + vsp * alpha * van - vap * vsn * (va2 + 1)) ...
    * beta * D1 ^ 2 * Cint1) / 0.3e1 + van * ((D1 ^ 2 * Sigmasur * vsn) + (CintD * beta * alpha) / 0.3e1)...
    * (beta * va2 * vap - vs2 * vsp)) / van / (vap * va2 * (D1 - 1) * (D1 + 1) * beta + vsp * (D1 ^ 2 + vs2)) / vsn / beta;

x = beta*ASEI/vs2/vsp-beta*BS/vsp/r^2-alpha*beta*Cintfun(r,t,k)/3/r^2/vsn;

end

