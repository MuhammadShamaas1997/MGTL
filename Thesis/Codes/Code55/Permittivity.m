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

epsr=0.9999;
a0=1e-4;%0.1mm
c0=2.99792458e8;%Speed of Light (m/s)
f0=c0/a0;%3000GHz
t0=1/f0;%0.33e-12 (s)
mu0=4*pi*(1e-7);% (H/m)
eps0=8.854187817e-12;% (F/m)

loglog(f,eps);
grid('off')
xlabel('Frequency (Hz)')
ylabel('Relative Permeability \mu_r')
axis([7e7 1e10 0.1 1e4])
