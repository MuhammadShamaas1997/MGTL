clc;clear all;
f=fopen('Flux.txt');
l=fgetl(f);
i=1;
while ischar(l)
    %%disp(l);
    text{i}=l;
    data{i}=sscanf(text{i},'%f , %f , %f');
    A1(i)=data{i}(1);
    A2(i)=data{i}(2);
    A3(i)=data{i}(3);
    l=fgetl(f);
    i=i+1;
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
S0=(I0*I0)/(eps0*c0*a0*a0);%//Poynting Vector


hold on;

subplot(3,1,1)
%A=A2);%/max(abs(A2));%Flux_in
%B=abs(A3);%/max(abs(A2));%Flux_out

plot(A1*f0,-A2*S0);
%axis([1e7 1e11 0 2e9]);
ylabel('|Sin| (W/m^2)');
xlabel('frequency (Hz)');
grid('on')

subplot(3,1,2)
plot(A1*f0,-A3*S0);
%axis([1e7 1e11 0 3e7]);
ylabel('|Sout| (W/m^2)');
xlabel('frequency (Hz)');
grid('on')

subplot(3,1,3)
C=10*log(A2./A3);
plot(f0*A1,C);
%axis([1e7 1e11 39.2 39.7]);
ylabel('Insertion Loss (dB)');
xlabel('frequency (Hz)');
grid('on')

