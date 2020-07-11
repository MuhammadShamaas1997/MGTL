clc;clear all;
%hold on;
%%
f=fopen('FieldEvolutionIn.txt');
File=1;
l=fgetl(f);
in=1;
Prev=1600*(File-1);
while ischar(l)
%while (in <= (257*81))
%for kj=1:81
    %%disp(l);
    text{in}=l;
    data{in}=sscanf(text{in},'%f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f');
    Hx(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(3)+1i*data{in}(4);
    Hy(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(5)+1i*data{in}(6);
    Hz(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(7)+1i*data{in}(8);
    
    Bx(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(9)+1i*data{in}(10);
    By(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(11)+1i*data{in}(12);
    Bz(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(13)+1i*data{in}(14);
    
    Ex(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(15)+1i*data{in}(16);
    Ey(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(17)+1i*data{in}(18);
    Ez(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(19)+1i*data{in}(20);
   
    Dx(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(21)+1i*data{in}(22);
    Dy(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(23)+1i*data{in}(24);
    Dz(data{in}(1)-Prev,(((data{in}(2)+15)*4)+1))=data{in}(25)+1i*data{in}(26);
    
    l=fgetl(f);
    in=in+1;
end


for t=1:10:128
plot(imag(Bz(t,:)));hold on;
pause(0.1);
end

% hold on;
% plot(0:(2/400):2,(real(Bx(170,:))*B0));
% plot(0:(2/400):2,(imag(Bx(170,:))*B0),'r');
% grid('on');
% xlabel('Distance (mm)');
% ylabel('|Bz (Wb/m^2)|');
% legend('Real','Imaginary')

epsr=1;
a0=1e-1;%0.1mm
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
sigmaD0=(epsr*eps0*c0)/a0;

T=t0;
Fs=1/T;
L=256/2;
L=2^nextpow2(L);
f=Fs/2*linspace(0,1,L/2+1);
FHxi=(fft(Hx(1:L,60),L));
FHxo=(fft(Hx(1:L,61),L));

Gamma=log(FHxo./FHxi)/(-(1/120)*(030*a0));
Gamma(1)=0;

figure;
subplot(2,1,1)
hold on;
plot(t0*(1:L),real(Hz(1:L,60)));
xlabel('Time (s)');ylabel('Real Hz(t)')
subplot(2,1,2)
hold on;
plot(t0*(1:L),imag(Hz(1:L,60)));
xlabel('Time (s)');ylabel('Imaginary Hz(t)')

figure
subplot(2,1,1)
plot(f(1:L/2+1),abs(FHxi(1:L/2+1)))
xlabel('Frequency (s)');ylabel('|Hz (jw)|')
%axis([0 5e9 0 2e-5])
subplot(2,1,2)
plot(f(1:L/2+1),angle(FHxi(1:L/2+1))*(180/pi))
xlabel('Frequency (s)');ylabel('\theta Hz (jw)')
%axis([0 5e9 -200 200])


figure;
subplot(3,1,1)
hold on;
plot(f(1:L/2+1),abs(real(Gamma(1:L/2+1))));
%hold on; plot(f(1:L/2+1),2*pi*f(1:L/2+1)/(3e8),'r');
xlabel('Frequency (Hz)')
ylabel('\alpha (Np.m^-^1)');
%axis([0 5e9 0 3e2])

subplot(3,1,2)
plot(f(1:L/2+1),abs(imag(Gamma(1:L/2+1))));
%hold on; plot(f(1:L/2+1),2*pi*f(1:L/2+1)/(3e8),'r');
ylabel('\beta (rad.m^-^1)');
xlabel('Frequency (Hz)')
%axis([0 5e9 0 5e2])

f=f';
subplot(3,1,3)
plot(f(1:L/2+1),abs(2*pi*f(1:L/2+1)./(imag(Gamma(1:L/2+1)))));
%hold on; plot(f(1:L/2+1),(3e8),'r');
ylabel('vp (m.s^-^1)');
xlabel('Frequency (Hz)')
%axis([0 5e9 0 2e8])

% X = 1/(4*sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)));

T=t0;
Fs=1/T;
L=256/2;
NFFT=2^nextpow2(L);
f=Fs/2*linspace(0,1,NFFT/2+1);
FEx=fft(Ex(1:L,60)*E0,NFFT)/L;
FHx=fft(Hx(1:L,60)*H0,NFFT)/L;
Z=FEx./FHx;

figure;
subplot(2,1,1)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
plot(f(1:NFFT/2+1),abs(Z(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|\eta| (Ohm)');
%axis([0 5e9 0 5000])

subplot(2,1,2)
plot(f(1:NFFT/2+1),abs(angle(Z(1:NFFT/2+1))*(180/pi)));
ylabel('\Theta \eta (Degree)');
xlabel('Frequency (Hz)');
%axis([0 5e9 0 200])



Zimag=abs(imag(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
Zreal=abs(real(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
Zabs=abs((Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
Zangle=abs(angle((Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)))*(180/pi));
Gm=abs(real(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
Rm=abs((Gm')*(-1).*(2*pi*f));%Reluctance
XCm=abs(imag(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
Yabs=abs((Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
Yangle=abs(angle((Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)))*(180/pi));

figure;
subplot(2,1,1)
plot(f,(Zabs(1:NFFT/2+1)));
ylabel('|Z| (Ohm/m)');
xlabel('Frequency (Hz)');
%axis([0 5e9 0 10e5])

subplot(2,1,2);
plot(f,(Zangle(1:NFFT/2+1)),'r');
ylabel('\theta Z (Degree)');
xlabel('Frequency (Hz)');
%axis([0 5e9 0 200])

% subplot(3,1,1);
% plot(f(1:NFFT/2+1),(Gm(1:NFFT/2+1)));
% ylabel('Conductance GL (S/m)');
% xlabel('Frequency (Hz)');
% % axis([0 5e11 -8.5 -7.5])
% 
% subplot(3,1,2);
% plot(f(1:NFFT/2+1),(Rm(1:NFFT/2+1)));
% ylabel('Reluctance Rmskin (1/H.m)');
% xlabel('Frequency (Hz)');
% % axis([0 5e11 0 3e13])
% 
% subplot(3,1,3);
% plot(f(1:NFFT/2+1),(XCm(1:NFFT/2+1)));
% ylabel('Susceptance XCL (S/m)');
% xlabel('Frequency (Hz)');
% % axis([0 5e11 0 1])

figure;
subplot(2,1,1)
plot(f,(Yabs(1:NFFT/2+1)));
ylabel('Y (1/Ohm.m)');
xlabel('Frequency (Hz)');
%axis([0 5e9 0 2])

subplot(2,1,2);
plot(f,(Yangle(1:NFFT/2+1)),'r');
ylabel('\theta Y (Degree)');
xlabel('Frequency (Hz)');
%axis([0 5e9 0 200])
