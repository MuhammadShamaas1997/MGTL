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

epsr=10;
a0=1e-4;%0.1mm
c0=2.99792458e8;%Speed of Light (m/s)
f0=c0/a0;%3000GHz
t0=1/f0;%0.33e-12 (s)
mu0=4*pi*(1e-7);% (H/m)
eps0=8.854187817e-12;% (F/m)
sigmaD0=(epsr*eps0*c0)/a0;

sigma=1e4;
wn=0.01;
gamma=2*pi;
muinf=1;
% chi2=chi3=1;
% conductivity=0;
% lorentzian(0.01,1,true);
ind=1;
for w=((1e8)/f0):((1e9)/f0):((1e12)/f0)
    mu(ind)=muinf+((sigma*wn*wn)/(wn*wn-w*w-1i*w*gamma));
    ind=ind+1;
end
w=((1e8)/f0):((1e9)/f0):((1e12)/f0);
w=w*f0;
loglog(w,real(mu))
hold on;
loglog(w,imag(mu),'r')
legend('Real','Imaginary');
%figure;
loglog(f,eps);
% hold on;
figure;
subplot(2,2,1);
g=(sqrt((1i*2*pi*w.*mu*mu0).*((5e-3)+1i*2*pi*w*10*eps0)));
plot(w,abs(real(g)));
subplot(2,2,1);
xlabel('Frequency (Hz)');ylabel('Attenuation Constant \alpha');
subplot(2,2,2);
plot(w,abs(imag(g)));
xlabel('Frequency (Hz)');ylabel('Phase Constant \beta');


eta=(sqrt((1i*2*pi*w.*mu*mu0)./((5e-3)+1i*2*pi*w*10*eps0)));
subplot(2,2,3);
plot(w,abs(eta));
xlabel('Frequency (Hz)');ylabel('Intrinsic Impedance \eta');
subplot(2,2,4);
plot(w,angle(eta)*(180/pi));
xlabel('Frequency (Hz)');ylabel('Impedance Angle \theta \eta');

XL=(abs(g.*eta));
XLm=(imag(g.*eta));
R=(real(g.*eta));
Gm=(real(g./eta));
%Rm=(Gm')*(-1).*(2*pi*w);%Reluctance
XCm=(imag(g./eta));

figure;
ld0=50;
ld=50;
ld=ld+1./(Gm+1i*XCm);
ld=(ld.*(1i*XLm))./(ld+(1i*XLm));
ld=ld+1./(Gm+1i*XCm);
plot(w,abs(ld0./ld))