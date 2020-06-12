clc;
clear all;
f=1e9;
mu0=4*pi*(1e-7);% (H/m)
eps0=8.854187817e-12;% (F/m)
indmu=1;
indsig=1;
for mur=1:(1e2):(1e4)
    indsig=1;
    for resistivity=(1e5):(1e9):(1e11)
        sigma=1/resistivity;
        Z(indmu,indsig)=sqrt((i*2*pi*f*mur*mu0)/(sigma+i*2*pi*f*eps0));
        
        indsig=indsig+1;
    end
    indmu=indmu+1;
end
hold on;
surf(1:(1e2):(1e4),(1e5):(1e9):(1e11),abs(Z))
grid on;
xlabel('Relative Permeability \mu_r')
ylabel('Electrical Resistivity \rho')
zlabel('Intrinsic Impedance \Omega')







% clc;
% clear all;
% sigma=1e9;
% mu0=4*pi*(1e-7);% (H/m)
% eps0=8.854187817e-12;% (F/m)
% indmu=1;
% indf=1;
% for mur=1:(1e2):(1e4)
%     indf=1;
%     for f=(1e5):(1e9):(1e11)
%         %sigma=1/resistivity;
%         Z(indmu,indf)=sqrt((i*2*pi*f*mur*mu0)/(sigma+i*2*pi*f*eps0));
%         
%         indf=indf+1;
%     end
%     indmu=indmu+1;
% end
% hold on;
% surf(1:(1e2):(1e4),(1e5):(1e9):(1e11),abs(Z))
% grid on;
% xlabel('Relative Permeability \mu_r')
% ylabel('Frequency (Hz)')
% zlabel('Intrinsic Impedance \Omega')