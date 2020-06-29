clc;clear all;
%hold on;
%%
f=fopen('F1.txt');
File=1;
l=fgetl(f);
in=1;
Prev=1600*(File-1);
while ischar(l)
%for kj=1:81
    %%disp(l);
    text{in}=l;
    data{in}=sscanf(text{in},'%f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f');
    Hx(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(3)+1i*data{in}(4);
    Hy(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(5)+1i*data{in}(6);
    Hz(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(7)+1i*data{in}(8);
    
    Bx(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(9)+1i*data{in}(10);
    By(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(11)+1i*data{in}(12);
    Bz(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(13)+1i*data{in}(14);
    
    Ex(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(15)+1i*data{in}(16);
    Ey(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(17)+1i*data{in}(18);
    Ez(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(19)+1i*data{in}(20);
   
    Dx(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(21)+1i*data{in}(22);
    Dy(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(23)+1i*data{in}(24);
    Dz(data{in}(1)-Prev,(((data{in}(2)+50)*4)+1))=data{in}(25)+1i*data{in}(26);
    
    l=fgetl(f);
    in=in+1;
end

% hold on;
% for ti=1:10:256
%     plot((abs(Bz(ti,:))));
%     pause(1)
% end

% hold on;
% plot(0:(2/400):2,(real(Bx(170,:))*B0));
% plot(0:(2/400):2,(imag(Bx(170,:))*B0),'r');
% grid('on');
% xlabel('Distance (mm)');
% ylabel('|Bz (Wb/m^2)|');
% legend('Real','Imaginary')

epsr=10;
a0=1e-2;%0.1mm
c0=2.99792458e8;%Speed of Light (m/s)
f0=c0/a0;%3000GHz
t0=1/f0;%0.33e-12 (s)
mu0=4*pi*(1e-7);% (H/m)
eps0=8.854187817e-12;% (F/m)
I0=1; %(A)
E0=I0/(a0*eps0*c0);%Electric Field
D0=I0/(a0*c0);%Electric Displacement Field
B0=I0/(a0*eps0*c0*c0);%Magnetic Field
H0=I0/(a0);%Magnetizing Field
sigmaD0=(epsr*eps0*c0)/a0;

hold on;

hold on;
T=t0;
Fs=1/T;
L=256;
L=2^nextpow2(L);
hold on

FHxi=(fft(Hx(1:L,2),L));
FHxi=FHxi(1:L/2+1);
f=Fs*(0:L/2+1);
f=f';

FHxo=(fft(Hx(1:L,3),L));
FHxo=FHxo(1:L/2+1);

Gamma=log(FHxo./FHxi)/(-(1/400)*(100*a0));

subplot(3,1,1)
hold on;
plot(f(2:L/2+1),real(Gamma(2:L/2+1)));
xlabel('Frequency (Hz)')
ylabel('\alpha (Np.m^-^1)');
axis([1e9 1e10 710 730])

subplot(3,1,2)
plot(f(2:L/2+1),abs(imag(Gamma(2:L/2+1))));
ylabel('\beta (rad.m^-^1)');
xlabel('Frequency (Hz)')
axis([1e9 1e10 215 220])

subplot(3,1,3)
plot(f(2:L/2+1),abs(2*pi*f(2:L/2+1)./(imag(Gamma(2:L/2+1)))));
ylabel('vp (m.s^-^1)');
xlabel('Frequency (Hz)')
axis([1e9 1e10 0 5e8])
% X = 1/(4*sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)));

T=t0;
Fs=1/T;
L=256;
NFFT=2^nextpow2(L);
FEx=fft(Ex(1:L,2)*E0,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
FHx=fft(Hx(1:L,2)*H0,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
Z=FEx./FHx;
%Z=Z*(377/16.5);

subplot(2,1,1)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
plot(f(1:NFFT/2+1),abs(Z(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Z_w (Ohm)|');
axis([1e9 1e10 114.5  114.7])

subplot(2,1,2)
plot(f(1:NFFT/2+1),angle(Z(1:NFFT/2+1))*(180/pi));
ylabel('\Theta Z_w (deg)');
xlabel('Frequency (Hz)');
axis([1e9 1e10 -6 -5.5])



XLm=(imag(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
R=(real(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
Gm=(real(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
Rm=(Gm')*(-1).*(2*pi*f);%Reluctance
XCm=(imag(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));

subplot(2,1,1);
plot(f,abs(XLm(1:NFFT/2+1)));
ylabel('XLm(H/m)');
xlabel('Frequency (Hz)');
axis([1e9 1e10 3.29e4 3.32e4])

subplot(2,1,2);
plot(f,abs(R(1:NFFT/2+1)));
ylabel('R');
xlabel('Frequency (Hz)');
%axis([0 1e11 -5e9 5e9])


subplot(3,1,1);
plot(f(1:NFFT/2+1),abs(Gm(1:NFFT/2+1)));
ylabel('Conductance GL (S/m)');
xlabel('Frequency (Hz)');
axis([1e9 1e10 6.35 6.5])

subplot(3,1,2);
plot(f(1:NFFT/2+1),abs(Rm(1:NFFT/2+1)));
ylabel('Reluctance Rmskin (1/H.m)');
xlabel('Frequency (Hz)');
axis([1e9 1e10 0 5e11])

subplot(3,1,3);
plot(f(1:NFFT/2+1),abs(XCm(1:NFFT/2+1)));
ylabel('Susceptance XCL (S/m)');
xlabel('Frequency (Hz)');
axis([1e9 1e10 1.24 1.275])
