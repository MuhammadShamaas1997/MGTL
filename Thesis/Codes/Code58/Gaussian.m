clc;clear all;

epsr=1;
a0=1e-4;%0.1mm
c0=3e8;%2.99792458e8;%Speed of Light (m/s)
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


w=1e-5;
t0=5;
t=0:0.001:10;
g=exp(-1i*w*t-((t-t0).*(t-t0))./(2*w*w));
% hold on;
% plot(t,abs(g));
% plot(t,real(g),'r')
% plot(t,imag(g),'c')

f=fopen('SourceFFT.txt');
l=fgetl(f);
in=1;
while ischar(l)
    text{in}=l;
    data{in}=sscanf(text{in},'%f %f %f');
    fr(in)=data{in}(1);
    FFT(in)=data{in}(2)+1i*data{in}(3);
    
    l=fgetl(f);
    in=in+1;
end
hold on;
plot(fr*f0,abs(FFT))
plot(fr*f0,abs(FFT),'r')
