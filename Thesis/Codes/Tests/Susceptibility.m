clc;
clear all;

eV_um_scale = 1.0/1.23984193;
Cu_plasma_frq = 10.83*eV_um_scale;
Cu_f0 = 0.575;
Cu_frq0 = 1e-10;
Cu_gam0 = 0.030*eV_um_scale;
Cu_sig0 = Cu_f0*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq0*Cu_frq0);
Cu_f1 = 0.061;
Cu_frq1 = 0.291*eV_um_scale;
Cu_gam1 = 0.378*eV_um_scale;
Cu_sig1 = Cu_f1*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq1*Cu_frq1);
Cu_f2 = 0.104;
Cu_frq2 = 2.957*eV_um_scale;      
Cu_gam2 = 1.056*eV_um_scale;
Cu_sig2 = Cu_f2*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq2*Cu_frq2);
Cu_f3 = 0.723;
Cu_frq3 = 5.300*eV_um_scale;      
Cu_gam3 = 3.213*eV_um_scale;
Cu_sig3 = Cu_f3*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq3*Cu_frq3);
Cu_f4 = 0.638;
Cu_frq4 = 11.18*eV_um_scale;
Cu_gam4 = 4.305*eV_um_scale;
Cu_sig4 = Cu_f4*(Cu_plasma_frq*Cu_plasma_frq)/(Cu_frq4*Cu_frq4);
sigma_Cu=Cu_sig0;

um_scale=1000;
NiFe_frq = 1/(0.0838297450980392*um_scale);
NiFe_gam = 1/(0.259381156903766*um_scale);
NiFe_sig = 1;
mu0=1000;
sigmab=0;
sigma0=NiFe_sig;
omega0=2*pi*NiFe_frq;
gamma0=NiFe_gam;
n=1;
for omega=0:0.0001:0.1
    mu(n)=(1+(i*sigmab/omega))*(mu0+(sigma0*omega0)/(omega0*omega0-omega*omega-i*omega*gamma0));
    n=n+1;
end
plot(0:0.0001:0.1,mu)


% eps0=1;
% sigmad=1000;
% n=1;
% sum=0;
% for omega=0:0.0001:0.01
%     sum=0;
%     sum=sum+(Cu_sig0*2*pi*Cu_frq0)/(2*pi*Cu_frq0*2*pi*Cu_frq0-omega*omega-i*omega*Cu_gam0);
%     sum=sum+(Cu_sig1*2*pi*Cu_frq1)/(2*pi*Cu_frq1*2*pi*Cu_frq1-omega*omega-i*omega*Cu_gam1);
%     sum=sum+(Cu_sig2*2*pi*Cu_frq2)/(2*pi*Cu_frq2*2*pi*Cu_frq2-omega*omega-i*omega*Cu_gam2);
%     sum=sum+(Cu_sig3*2*pi*Cu_frq3)/(2*pi*Cu_frq3*2*pi*Cu_frq3-omega*omega-i*omega*Cu_gam3);
%     sum=sum+(Cu_sig4*2*pi*Cu_frq4)/(2*pi*Cu_frq4*2*pi*Cu_frq4-omega*omega-i*omega*Cu_gam4);
%     
%     eps(n)=(1+(i*sigmad/omega))*(eps0+sum);
%     n=n+1;
% end
% plot(0:0.0001:0.01,abs(eps))
