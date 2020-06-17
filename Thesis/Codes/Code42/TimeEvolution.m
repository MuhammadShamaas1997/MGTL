clc;clear all;
f=fopen('TimeEvolution.txt');
l=fgetl(f);
i=1;
while ischar(l)
    %%disp(l);
    text{i}=l;
    data{i}=sscanf(text{i},'%f , %f , %f , %f , %f , %f , %f , %f');
    A1(i)=data{i}(1);
    A2(i)=data{i}(2);
    A3(i)=data{i}(3);
    A4(i)=data{i}(4);
    A5(i)=data{i}(5);
    A6(i)=data{i}(6);
    A7(i)=data{i}(7);
    A8(i)=data{i}(8);
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
V0=I0/(eps0*c0);%Electric Field

hold on;
A1=A1*V0;A2=A2*V0;
A3=A3*I0;A4=A4*I0;
A5=A5*I0;A6=A6*I0;
A7=A7*V0;A8=A8*V0;

t=100*t0*(1:length(A1));
subplot(2,1,1)
hold on;plot(t,A1);plot(t,A2,'r');%Im
ylabel('Magnetic Current (V)');
xlabel('Time');
legend('Re','Im');

subplot(2,1,2)
hold on;plot(t,A3);plot(t,A4,'r');%Vm
ylabel('Magnetic Voltage (A)');
xlabel('Time');
legend('Re','Im');

% subplot(2,2,3)
% hold on;plot(t,A5);plot(t,A6,'r');%Ie
% ylabel('Electric Current (A)');
% xlabel('Time')
% legend('Re','Im');
% 
% subplot(2,2,4)
% hold on;plot(t,A7);plot(t,A8,'r');%Ve
% ylabel('Electric Voltage (V)');
% xlabel('Time')
% legend('Re','Im');
% 
