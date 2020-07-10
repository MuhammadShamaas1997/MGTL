clc;clear all;
f=fopen('FieldEvolutionIn.txt');
l=fgetl(f);
in=1;
while ischar(l)
    %%disp(l);
    text{in}=l;
    data{in}=sscanf(text{in},'%f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f');
    A1(in)=data{in}(1);
    A2(in)=data{in}(2);
    A3(in)=data{in}(3);
    A4(in)=data{in}(4);
    A5(in)=data{in}(5);
    A6(in)=data{in}(6);
    A7(in)=data{in}(7);
    A8(in)=data{in}(8);
    A9(in)=data{in}(9);
    A10(in)=data{in}(10);
    A11(in)=data{in}(11);
    A12(in)=data{in}(12);
    A13(in)=data{in}(13);
    A14(in)=data{in}(14);
    A15(in)=data{in}(15);
    A16(in)=data{in}(16);
    A17(in)=data{in}(17);
    A18(in)=data{in}(18);
    A19(in)=data{in}(19);
    A20(in)=data{in}(20);
    A21(in)=data{in}(21);
    A22(in)=data{in}(22);
    A23(in)=data{in}(23);
    A24(in)=data{in}(24);
    l=fgetl(f);
    in=in+1;
end

f=fopen('FieldEvolutionOut.txt');
l=fgetl(f);
in=1;
while ischar(l)
    %%disp(l);
    text{in}=l;
    data{in}=sscanf(text{in},'%f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f');
    A1o(in)=data{in}(1);
    A2o(in)=data{in}(2);
    A3o(in)=data{in}(3);
    A4o(in)=data{in}(4);
    A5o(in)=data{in}(5);
    A6o(in)=data{in}(6);
    A7o(in)=data{in}(7);
    A8o(in)=data{in}(8);
    A9o(in)=data{in}(9);
    A10o(in)=data{in}(10);
    A11o(in)=data{in}(11);
    A12o(in)=data{in}(12);
    A13o(in)=data{in}(13);
    A14o(in)=data{in}(14);
    A15o(in)=data{in}(15);
    A16o(in)=data{in}(16);
    A17o(in)=data{in}(17);
    A18o(in)=data{in}(18);
    A19o(in)=data{in}(19);
    A20o(in)=data{in}(20);
    A21o(in)=data{in}(21);
    A22o(in)=data{in}(22);
    A23o(in)=data{in}(23);
    A24o(in)=data{in}(24);
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


hold on;

subplot(4,1,1)
hold on;%plot(mag3(mag(A1,A2),mag(A3,A4),mag(A5,A6)));%H
for m=1:length(A5)
        Hyi(m)=A3(m)+i*A4(m);
        Hyi(m)=Hyi(m)*H0;
end
T=100*t0;
Fs=1/T;
L=length(A1);
NFFT=2^nextpow2(L);
FHyi=fft(Hyi,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
%plot(f,2*abs(FHyi(1:NFFT/2+1)));
% plot(f,2*angle(YHz(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Hyi(f)|');
axis([2e9 15e9 0  1e-3])
%axis([])

%subplot(2,2,2)
%hold on;%plot(mag3(mag(A7,A8),mag(A9,A10),mag(A11,A12)));%B

subplot(4,1,2)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
for m=1:length(A3o)
        Hyo(m)=A3o(m)+i*A4o(m);
        Hyo(m)=Hyo(m)*H0;
end
T=100*t0;
Fs=1/T;
L=length(A1);
NFFT=2^nextpow2(L);
FHyo=fft(Hyo,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);

%plot(f,2*abs(FHyo(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Hyo(A/m)|');
axis([2e9 15e9 0  1e-3])
Gamma=log(FHyo./FHyi)/(-0.5e-3);

subplot(2,1,1)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
loglog(f,2*real(Gamma(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('\alpha (m^-^1)');
%axis([0.1e9 15e9 0  10000])

% subplot(2,2,4)
% hold on;plot(mag3(mag(A19,A20),mag(A21,A22),mag(A23,A24)));%D

subplot(2,1,2)
plot(f,2*imag(Gamma(1:NFFT/2+1))*(180/pi),'.-');
ylabel('\beta (m^-^1)');
xlabel('Frequency (Hz)')
%axis([0.1e9 15e9 0 10e5])