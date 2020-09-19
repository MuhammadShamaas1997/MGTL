clc;clear all;

fi=fopen('Permeability.txt');
l=fgetl(fi);
in=1;
while ischar(l)
%for kj=1:81
    %%disp(l);
    text{in}=l;
    data{in}=sscanf(text{in},'%f %f');
    f(in)=(data{in}(1));
    eps(in)=data{in}(2);
    
    l=fgetl(fi);
    in=in+1;
end

epsr=0.9999;
a0=1e-2;%0.1mm
c0=2.99792458e8;%Speed of Light (m/s)
f0=c0/a0;%3000GHz
t0=1/f0;%0.33e-12 (s)
mu0=4*pi*(1e-7);% (H/m)
eps0=8.854187817e-12;% (F/m)
I0=1;
E0=I0/(a0*eps0*c0);
D0=I0/(a0*c0);
B0=I0/(a0*eps0*c0*c0);
H0=I0/(a0);
sigmaD0=(epsr*eps0*c0)/a0;
J0=I0/(a0*a0);
u0=(I0*I0)/(eps0*c0*c0*a0*a0);
S0=(I0*I0)/(eps0*c0*a0*a0);
Sc0=1/(c0);
sig0=-20;


semilogx(f/(a0/(1e-2)),eps);xlim([10e5 1e7]);
grid('off')
xlabel('Frequency (Hz)')
ylabel('Relative Permeability \mu_r')
title('MEEP')
%axis([1e8 1e9 0 1e4])

% muinf=1;gamma=.01/(8*4*pi*pi*pi*pi*pi);fn=0.01/(4*pi*pi);
% f=0:1e-7:(1/3);
% fi=f*(4*pi*pi);
% sigma=-100*(4*pi*pi);
% kwi=(fn.*fn-fi.*fi-1i.*gamma.*fi);
% b0=1.0;
% kwi2=(fi.*fi.*b0.*b0)./kwi;
% mur=muinf+(sigma.*fn.*fn)./(kwi-kwi2);

fn=(2/3)*(1e-5);
sigma=-1e4;
muinf=1;gamma=(.1e-5)/(2*pi);
f=0:1e-7:(1/30);
fi=f;
kwi=(fn.*fn-fi.*fi-1i.*gamma.*fi);
b0=1e-5;
kwi2=(fi.*fi.*b0.*b0)./kwi;
mur=muinf+(sigma.*fn.*fn)./(kwi-kwi2);

% figure;
% subplot(2,1,1);semilogx(f*f0,real(mur));title('Formula');xlim([10e5 1e7]);
% subplot(2,1,2);semilogx(f*f0,imag(mur));xlim([10e5 1e7]);

f=1e3:1e3:1e9;
nomf=((f)./(0.2e6));
cf=(1+nomf.*nomf);
mur2=mu0+((10000*mu0)./cf)-(1i.*nomf.*(10000*mu0))./cf;


figure;
subplot(2,1,1);hold on;
f=0:1e-7:(1/30);
semilogx(f*f0,real(mur));title('Formula');xlim([10e5 1e7]);
f=1e3:1e3:1e9;
semilogx(f,real(mur2)/mu0,'r');title('Literature');xlim([10e5 1e7]);
subplot(2,1,2);hold on;
f=0:1e-7:(1/30);
semilogx(f*f0,imag(mur));xlim([10e5 1e7]);
f=1e3:1e3:1e9;
semilogx(f,imag(mur2)/mu0,'r');xlim([10e5 1e7]);

