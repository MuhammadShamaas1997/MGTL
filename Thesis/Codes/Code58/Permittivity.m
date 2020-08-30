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

semilogx(f/(a0/1e-4),eps);
grid('off')
xlabel('Frequency (Hz)')
ylabel('Relative Permeability \mu_r')
title('MEEP')
%axis([1e8 1e9 0 1e4])

muinf=1;gamma=.01/(8*4*pi*pi*pi*pi*pi);fn=0.01/(4*pi*pi);
sigma=-20*(4*pi*pi);
f=0:1e-7:(1/3);
fi=f*(4*pi*pi);
%mur=muinf+(sigma.*fn.*fn)./(-fi.*fi-1i.*gamma.*f);
sigma=-100*(4*pi*pi);
mur=muinf+(sigma.*fn.*fn)./(fn.*fn-fi.*fi-1i.*gamma.*fi);
figure;
subplot(2,1,1);semilogx(f*f0,real(mur));title('Formula');
subplot(2,1,2);semilogx(f*f0,imag(mur));

f=1e3:1e3:1e9;
nomf=((f)./(0.2e6));
cf=(1+nomf.*nomf);
mur2=mu0+((10000*mu0)./cf)-(1i.*nomf.*(10000*mu0))./cf;


figure;
subplot(2,1,1);semilogx(f,real(mur2)/mu0);title('Literature');
subplot(2,1,2);semilogx(f,imag(mur2)/mu0);


% f=f*f0;mur=mur*mu0;sigma=5e-3;eps=(1*eps0)-1i*(sigma./(2*pi*f));
% gamma=1i.*2.*pi.*f.*sqrt(mur.*eps);
% eta=sqrt(mur./eps);
% Z=gamma.*eta;Y=gamma./eta;
% alpha=real(gamma);beta=imag(gamma);vp=(2.*pi.*f)./beta;
% 
% figure;hold on;
% subplot(2,1,1);semilogx(f,abs(Z));title('|Z|');
% subplot(2,1,2);semilogx(f,angle(Z)*(180/pi));title('theta Z');
% 
% figure;hold on;
% subplot(2,1,1);semilogx(f,abs(Y));title('|Y|');
% subplot(2,1,2);semilogx(f,angle(Y)*(180/pi));title('theta Y');
% 
% figure;hold on;
% subplot(2,1,1);semilogx(f,alpha);title('alpha');
% subplot(2,1,2);semilogx(f,beta);title('beta');
% 
% figure;hold on;
% subplot(2,1,1);semilogx(f,abs(eta));title('\eta');
% subplot(2,1,2);semilogx(f,angle(eta)*(180/pi));title('theta \eta');
% 
% figure;
% semilogx(f,vp);title('vp');
