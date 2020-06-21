clc;clear all;
f=fopen('Skin.txt');
l=fgetl(f);
i=1;
while ischar(l)
    %%disp(l);
    text{i}=l;
    data{i}=sscanf(text{i},'%f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f');
    A1(i)=data{i}(1);
    A2(i)=data{i}(2);
    A3(i)=data{i}(3);
    A4(i)=data{i}(4);
    A5(i)=data{i}(5);
    A6(i)=data{i}(6);
    A7(i)=data{i}(7);
    A8(i)=data{i}(8);
    A9(i)=data{i}(9);
    A10(i)=data{i}(10);
    A11(i)=data{i}(11);
    A12(i)=data{i}(12);
    A13(i)=data{i}(13);
    A14(i)=data{i}(14);
    A15(i)=data{i}(15);
    A16(i)=data{i}(16);
    A17(i)=data{i}(17);
    A18(i)=data{i}(18);
    A19(i)=data{i}(19);
    A20(i)=data{i}(20);
    A21(i)=data{i}(21);
    A22(i)=data{i}(22);
    A23(i)=data{i}(23);
    A24(i)=data{i}(24);
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

hold on;

subplot(2,2,1)
hold on;plot(-5:0.001:5-0.001,log(H0*mag3(mag(A1,A2),mag(A3,A4),mag(A5,A6))));%Hx
ylabel('log |H (A/m)|');
xlabel('Length (mm)');

subplot(2,2,2)
hold on;plot(-5:0.001:5-0.001,log(B0*mag3(mag(A7,A8),mag(A9,A10),mag(A11,A12))));%Hx
ylabel('log |B (Wb/m^2)|');
xlabel('Length (mm)');

subplot(2,2,3)
hold on;plot(-5:0.001:5-0.001,log(E0*mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18))));%Hx
ylabel('log |E (V/m)|');
xlabel('Length (mm)');

subplot(2,2,4)
hold on;plot(-5:0.001:5-0.001,log(D0*mag3(mag(A19,A20),mag(A21,A22),mag(A23,A24))));%Hx
ylabel('log |D (C/m^2)|');
xlabel('Length (mm)');


