clc;clear all
epsr=1.0;

%SI Conversion Factors
a0=1e-3;%1mm
c0=2.99792458e8;%Speed of Light (m/s)
f0=c0/a0;%300GHz
t0=1/f0;%0.33e-11 (s)
mu0=4*pi*(1e-7);% (H/m)
eps0=8.854187817e-12;% (F/m)
I0=1;%(A)
E0=I0/(a0*eps0*c0);%Electric Field
D0=I0/(a0*c0);%Electric Displacement Field
B0=I0/(a0*eps0*c0*c0);%Magnetic Field
H0=I0/(a0);%Magnetizing Field
sigmaD0=(epsr*eps0*c0)/a0;%Electric Conductivity
J0=I0/(a0*a0);%Electric Current Density
u0=(I0*I0)/(eps0*c0*c0*a0*a0);%Energy Density
S0=(I0*I0)/(eps0*c0*a0*a0);%Poynting Vector
Sc0=1/(c0);%Courant Factor
