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
sigmaD0=(epsr*eps0*c0)/a0;%Electric Conductivity

hold on;

subplot(3,1,1)
% A=abs(A2)/max(abs(A2));%Flux_in
% B=abs(A2)/max(abs(A3));%Flux_out

subplot(3,1,1)
plot(A1*f0,abs(-A2*S0));
axis([0 4e9 S0*0.01 S0*40]);
ylabel('|Sin| (W/m^2)');
xlabel('frequency (Hz)');


subplot(3,1,2)
plot(A1*f0,abs(-A3*S0));
axis([0 4e9 S0*0.01 S0*25]);
ylabel('|Sout| (W/m^2)');
xlabel('frequency (Hz)');


subplot(3,1,3)
C=log(A2./A3);
C=C(1:5:1000);
A1=A1(1:5:1000);
plot((f0*A1),(abs(C)));
axis([0 4e9 -1 5]);
ylabel('Insertion Loss');
xlabel('frequency (Hz)');
