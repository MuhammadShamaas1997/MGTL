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
    Hx(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(3)+1i*data{in}(4);
    Hy(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(5)+1i*data{in}(6);
    Hz(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(7)+1i*data{in}(8);
    
    Bx(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(9)+1i*data{in}(10);
    By(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(11)+1i*data{in}(12);
    Bz(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(13)+1i*data{in}(14);
    
    Ex(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(15)+1i*data{in}(16);
    Ey(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(17)+1i*data{in}(18);
    Ez(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(19)+1i*data{in}(20);
    
    Dx(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(21)+1i*data{in}(22);
    Dy(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(23)+1i*data{in}(24);
    Dz(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(25)+1i*data{in}(26);
    
    l=fgetl(f);
    in=in+1;
end

epsr=10;
a0=1e-4;%0.1mm
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

subplot(4,1,1)
hold on;
T=t0;
Fs=1/T;
L=45;
L=2^nextpow2(L);
% NFFT=2^nextpow2(L);
% FHxi=(fft(Hx(1:L,5)*H0,NFFT)/L)();
% FHxi(2:end-1)=
% f=Fs/2*linspace(0,1,NFFT/2+1);

FHxi=(fft(Hx(1:L,2)*H0,L));
FHxi=FHxi(1:L/2+1);
%FHxi(2:end-1)=2*FHxi(2:end-1)
f=Fs*(0:L/2+1);
f=f';

FHxo=(fft(Hx(1:L,12)*H0,L));
FHxo=FHxo(1:L/2+1);
%FHxo(2:end-1)=2*FHxo(2:end-1);

% hold on;
% T=t0;
% Fs=1/T;
% L=45;
% NFFT=2^nextpow2(L);
% FHxo=fft(Hx(1:L,10)*H0,NFFT)/L;
% f=Fs/2*linspace(0,1,NFFT/2+1);

%plot(f,2*abs(FHyo(1:NFFT/2+1)));
% xlabel('Frequency (Hz)')
% ylabel('|Hyo(A/m)|');
% axis([2e9 15e9 0  1e-3])

Gamma=log(FHxo./FHxi)/(-(10/80)*(2e-3));

subplot(3,1,1)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
plot(f(2:L/2+1),real(Gamma(2:L/2+1)));
xlabel('Frequency (Hz)')
ylabel('\alpha (Np.m^-^1)');
axis([0 5e11 1.32e5 1.38e5])

subplot(3,1,2)
plot(f(2:L/2+1),abs(imag(Gamma(2:L/2+1))));
ylabel('\beta (rad.m^-^1)');
xlabel('Frequency (Hz)')
axis([0 5e11 1e4 1.3e4])

subplot(3,1,3)
plot(f(2:L/2+1),f(2:L/2+1)./abs(imag(Gamma(2:L/2+1))));
ylabel('vp (m.s^-^1)');
xlabel('Frequency (Hz)')
axis([0 5e11 0 5e7])
% X = 1/(4*sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)));

T=t0;
Fs=1/T;
L=45;
NFFT=2^nextpow2(L);
FEx=fft(Ex(1:L,2)*E0,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
FHx=fft(Hx(1:L,2)*H0,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
Z=FEx./FHx;


subplot(2,1,1)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
plot(f(2:NFFT/2+1),2*abs(Z(2:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Z_w (Ohm)|');
axis([0 5e11 3.3e4  3.45e4])

subplot(2,1,2)
plot(f(2:NFFT/2+1),angle(Z(2:NFFT/2+1))*(180/pi));
ylabel('\Theta Z_w (deg)');
xlabel('Frequency (Hz)');
axis([0 5e11 170 190])


XLm=(imag(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
Gm=(real(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
Rm=(Gm')*(-1).*(2*pi*f);%Reluctance
XCm=(imag(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));

subplot(4,1,1);
plot(f,(XLm(1:NFFT/2+1)));
ylabel('XLm(H/m)');
xlabel('Frequency (Hz)');
%axis([0 1e11 -5e9 5e9])

subplot(3,1,1);
plot(f(2:NFFT/2+1),(Gm(2:NFFT/2+1)));
ylabel('Conductance GL (S/m)');
xlabel('Frequency (Hz)');
axis([0 5e11 -8.5 -7.5])

subplot(3,1,2);
plot(f(2:NFFT/2+1),(Rm(2:NFFT/2+1)));
ylabel('Reluctance Rmskin (1/H.m)');
xlabel('Frequency (Hz)');
axis([0 5e11 0 3e13])

subplot(3,1,3);
plot(f(2:NFFT/2+1),(XCm(2:NFFT/2+1)));
ylabel('Susceptance XCL (S/m)');
xlabel('Frequency (Hz)');
axis([0 5e11 0 1])
