% clc;
% clear all;
% f=1e9;
% mu0=4*pi*(1e-7);% (H/m)
% eps0=8.854187817e-12;% (F/m)
% indmu=1;
% indsig=1;
% for mur=1:(1e2):(1e4)
%     indsig=1;
%     for resistivity=(1e5):(1e9):(1e11)
%         sigma=1/resistivity;
%         Z(indmu,indsig)=sqrt((i*2*pi*f*mur*mu0)/(sigma+i*2*pi*f*eps0));
%         
%         indsig=indsig+1;
%     end
%     indmu=indmu+1;
% end
% hold on;
% surf(1:(1e2):(1e4),(1e5):(1e9):(1e11),abs(Z))
% grid on;
% xlabel('Relative Permeability \mu_r')
% ylabel('Electrical Resistivity \rho')
% zlabel('Intrinsic Impedance \Omega')







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










epsr=0.9999;
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


hold on;

subplot(4,1,1)
hold on;%plot(mag3(mag(A1,A2),mag(A3,A4),mag(A5,A6)));%H
for m=1:length(A1)
        Iz(m)=A1(m)+i*A2(m);
end
T=100*t0;
Fs=1/T;
L=length(A1);
NFFT=2^nextpow2(L);
FIz=fft(Iz,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
plot(f,2*abs(FIz(1:NFFT/2+1)));
% plot(f,2*angle(YHz(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Iz(f)|');
axis([0 5e9 0  3e4])


%subplot(2,2,2)
%hold on;%plot(mag3(mag(A7,A8),mag(A9,A10),mag(A11,A12)));%B

subplot(4,1,2)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
for m=1:length(A3)
        Vz(m)=A3(m)+i*A4(m);
end
T=100*t0;
Fs=1/T;
L=length(A1);
NFFT=2^nextpow2(L);
FVz=fft(Vz,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
plot(f,2*abs(FVz(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Ez(f)|');
axis([0 5e9 0 0.04])

Z=FIz./FVz;

subplot(4,1,3)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
plot(f,2*abs(Z(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Z(f)|');
axis([0 5e9 -1e8 1e8])

% subplot(2,2,4)
% hold on;plot(mag3(mag(A19,A20),mag(A21,A22),mag(A23,A24)));%D

subplot(4,1,4)
plot(f,2*imag(Z(1:NFFT/2+1))*(180/pi),'.-');
ylabel('\Theta Z(f)');
xlabel('Frequency (Hz)')
%axis([0 5e9 -200  200])