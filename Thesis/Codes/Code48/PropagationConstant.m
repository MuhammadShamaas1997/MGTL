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
% for m=1:length(A5)
%         Hyi(m)=A3(m)+i*A4(m);
%         Hyi(m)=Hyi(m)*H0;
% end
T=t0;
Fs=1/T;
L=45;
NFFT=2^nextpow2(L);
FHxi=fft(Hx(1:L,1)*H0,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
%plot(f,2*abs(FHyi(1:NFFT/2+1)));
% plot(f,2*angle(YHz(1:NFFT/2+1)));
% xlabel('Frequency (Hz)')
% ylabel('|Hyi(f)|');
% axis([2e9 15e9 0  1e-3])
%axis([])

%subplot(2,2,2)
%hold on;%plot(mag3(mag(A7,A8),mag(A9,A10),mag(A11,A12)));%B

% subplot(4,1,2)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
% for m=1:length(A3o)
%         Hyo(m)=A3o(m)+i*A4o(m);
%         Hyo(m)=Hyo(m)*H0;
% end
T=t0;
Fs=1/T;
L=45;
NFFT=2^nextpow2(L);
FHxo=fft(Hx(1:L,41)*H0,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);

%plot(f,2*abs(FHyo(1:NFFT/2+1)));
% xlabel('Frequency (Hz)')
% ylabel('|Hyo(A/m)|');
% axis([2e9 15e9 0  1e-3])

Gamma=log(FHxo./FHxi)/(-10e-3);

subplot(2,1,1)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
plot(f,abs(2*real(Gamma(1:NFFT/2+1))));
xlabel('Frequency (Hz)')
ylabel('\alpha (m^-^1)');
axis([0 1e10 1.157e4  1.16e4])

% subplot(2,2,4)
% hold on;plot(mag3(mag(A19,A20),mag(A21,A22),mag(A23,A24)));%D

subplot(2,1,2)
plot(f,abs(2*imag(Gamma(1:NFFT/2+1))),'.-');
ylabel('\beta (m^-^1)');
xlabel('Frequency (Hz)')
axis([0 1e10 550 650])






T=t0;
Fs=1/T;
L=45;
NFFT=2^nextpow2(L);
FEx=fft(Ex(1:L,1)*E0,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
FHx=fft(Hx(1:L,1)*H0,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
Z=FEx./FHx;
