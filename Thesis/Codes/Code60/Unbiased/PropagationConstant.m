clc;clear all;
%hold on;
%%
f=fopen('FieldEvolutionIn.txt');
File=1;
l=fgetl(f);
in=1;
Prev=1600*(File-1);
while ischar(l)
%while (in <= (518*300))
%for kj=1:81
    %%disp(l);
    text{in}=l;
    
    data{in}=sscanf(text{in},'%f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f');
    dis=round(((data{in}(2)+15)*4)+1);
    Hx(data{in}(1)-Prev,dis)=data{in}(3)+1i*data{in}(4);
    Hy(data{in}(1)-Prev,dis)=data{in}(5)+1i*data{in}(6);
    Hz(data{in}(1)-Prev,dis)=data{in}(7)+1i*data{in}(8);
    
    Bx(data{in}(1)-Prev,dis)=data{in}(9)+1i*data{in}(10);
    By(data{in}(1)-Prev,dis)=data{in}(11)+1i*data{in}(12);
    Bz(data{in}(1)-Prev,dis)=data{in}(13)+1i*data{in}(14);
    
    Ex(data{in}(1)-Prev,dis)=data{in}(15)+1i*data{in}(16);
    Ey(data{in}(1)-Prev,dis)=data{in}(17)+1i*data{in}(18);
    Ez(data{in}(1)-Prev,dis)=data{in}(19)+1i*data{in}(20);
   
    Dx(data{in}(1)-Prev,dis)=data{in}(21)+1i*data{in}(22);
    Dy(data{in}(1)-Prev,dis)=data{in}(23)+1i*data{in}(24);
    Dz(data{in}(1)-Prev,dis)=data{in}(25)+1i*data{in}(26);
    
    l=fgetl(f);
    in=in+1;
end

epsr=1;
a0=1e-2;%0.1mm
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


for ind=129:500
Bx(ind,:)=0;By(ind,:)=0;Bz(ind,:)=0;
Hx(ind,:)=0;Hy(ind,:)=0;Hz(ind,:)=0;
Dx(ind,:)=0;Dy(ind,:)=0;Dz(ind,:)=0;
Ex(ind,:)=0;Ey(ind,:)=0;Ez(ind,:)=0;
end

L2=32768*4;
%Bx(L2,121)=0;By(L2,121)=0;Bz(L2,121)=0;
Hx(L2,121)=0;Hy(L2,121)=0;Hz(L2,121)=0;
%Dx(L2,:)=0;Dy(L2,:)=0;Dz(L2,:)=0;
Ex(L2,121)=0;Ey(L2,121)=0;Ez(L2,121)=0;

% hold on;
% for stp=1:500
%     plot(real(Hx(stp,:)));
%     pause(.1);
% end

%Sx=Ex(1:100,:).*conj(Hx(1:100,:));
%Sy=Ey(1:100,:).*conj(Hy(1:100,:));
Sz=Ex(1:100,:).*conj(Hy(1:100,:))-Ey(1:100,:).*conj(Hx(1:100,:));
%S=abs(mag3(Sx,Sy,Sz));

figure;surf(abs(real(Sz(1:100,:))));hold on;
% figure;surf(abs(imag(Sz(1:200,:))));hold on;

figure;
hold on
for dist=1:121
plot3([0 real(Hx(90,dist))*H0],[dist dist+1], [0 0]);
plot3([0 0], [dist dist+1],[0 real(Hy(90,dist))*H0],'r');
%axis([-1e-4 1e-4 0 122 -1e-4 1e-4])
end
ylabel('Time Step n');
xlabel('Hx');
zlabel('Hy');
legend('Hx','Hy');
title('Real Field');


% hold on;
% plot(0:(2/400):2,(real(Bx(170,:))*B0));
% plot(0:(2/400):2,(imag(Bx(170,:))*B0),'r');
% grid('on');
% xlabel('Distance (mm)');
% ylabel('|Bz (Wb/m^2)|');
% legend('Real','Imaginary')

epsr=1;
a0=1e-2;%0.1mm
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
fmin=1e6;
fmax=1e9;

muinf=1;
gamma=(-0.33e-2)/(2*pi);
fn=(5.8e-5);
sigma=1e4;
fr=((fmin)/f0):1e-5:((fmax)/f0); 
fi=fr;
kwi=(fn.*fn-fi.*fi-1i.*gamma.*fi);
b0=0;
kwi2=(fi.*fi.*b0.*b0)./kwi;
mur=muinf+(sigma.*fn.*fn)./(kwi-kwi2);

fr2=fr*f0;
mur2=mur*mu0;
sigma=5e-3;
eps=(10*eps0)-1i*(sigma./(2*pi*fr2));


gammaor=1i.*2.*pi.*fr2.*sqrt(mur2.*eps);
etaor=sqrt(mur2./eps);
Zor=gammaor.*etaor;
Yor=gammaor./etaor;
alphaor=real(gammaor);
betaor=imag(gammaor);
vpor=(2.*pi.*fr2)./betaor;


figure;hold on;
subplot(2,1,1);semilogx(fr2,real(mur2)/mu0);ylabel('Real \mu_r');xlim([fmin fmax]);
title('Relative Permeability \mu_r');
subplot(2,1,2);semilogx(fr2,imag(mur2)/mu0);ylabel('Imaginary \mu_r');xlim([fmin fmax]);
xlabel('Frequency (Hz)')

figure
T=t0/4;
Fs=1/T;
L=L2;

dobs=1;
obs=90;
L=2^nextpow2(L);
f=Fs/2*linspace(0,1,L/2+1);

FHxi=(fft((Ey(1:L,obs)),L));
FHxo=(fft((Ey(1:L,obs+dobs)),L));
Gamma=(log(FHxo./FHxi))/(-(dobs/4)*(a0));
Gamma(1)=0;


figure;
subplot(2,1,1);
loglog(f(1:L/2+1),abs(real(Gamma(1:L/2+1))));
ylabel('\alpha (Np.m^-^1)');xlim([fmin fmax]);
hold on;loglog(fr2,alphaor,'r');ylabel('\alpha (Np.m^-^1)');xlim([fmin fmax]);
title('Propagation Constant \gamma');legend('Simulation','Theory');

subplot(2,1,2)
loglog(f(1:L/2+1),abs(imag(Gamma(1:L/2+1))));
ylabel('\beta (rad.m^-^1)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
hold on;loglog(fr2,betaor,'r');ylabel('\beta (rad.m^-^1)');xlim([fmin fmax]);
legend('Simulation','Theory');

figure
f=f';
loglog(f(1:L/2+1),abs(2*pi*f(1:L/2+1)./(imag(Gamma(1:L/2+1)))));
%hold on; plot(f(1:L/2+1),(3e8),'r');
ylabel('vp (m.s^-^1)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
hold on;loglog(fr2,vpor,'r');title('vp');xlim([fmin fmax]);


obs=25;
NFFT=2^nextpow2(L);
f=Fs/2*linspace(0,1,NFFT/2+1);
FEx=fft((Ey(1:L,obs)),NFFT);
FHx=fft((Hy(1:L,obs)),NFFT);
Z=FEx./FHx;
Z=Z*377;

figure;
subplot(2,1,1)
loglog(f(1:NFFT/2+1),abs(Z(1:NFFT/2+1)));
ylabel('|\eta| (Ohm)');xlim([fmin fmax]);
hold on;
loglog(fr2,abs(etaor),'r');title('Intrinsic Impedance \eta');xlim([fmin fmax]);
legend('Simulation','Theory');

subplot(2,1,2)
semilogx(f(1:NFFT/2+1),abs(angle(Z(1:NFFT/2+1))*(180/pi)));
ylabel('\Theta \eta (Degree)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
hold on;semilogx(fr2,abs(angle(etaor))*(180/pi),'r');xlim([fmin fmax]);ylim([-200 200]);
legend('Simulation','Theory');


Zimag=abs(imag(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
Zreal=abs(real(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
Zabs=abs((Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
Zangle=abs(angle((Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)))*(180/pi));
Gm=abs(real(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
Rm=abs((Gm')*(-1).*(2*pi*f));%Reluctance
XCm=abs(imag(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
Yabs=abs((Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
Yangle=(angle((Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)))*(180/pi));

figure;
subplot(2,1,1)
loglog(f,(Zabs(1:NFFT/2+1)));title('Transverse Impedance Z');
ylabel('|Z| (Ohm/m)');xlim([fmin fmax]);
hold on;
loglog(fr2,abs(Zor),'r');xlim([fmin fmax]);
legend('Simulation','Theory');

subplot(2,1,2);
semilogx(f,(Zangle(1:NFFT/2+1)));
ylabel('\theta Z (Degree)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
hold on;
semilogx(fr2,angle(Zor)*(180/pi),'r');xlim([fmin fmax]);ylim([-200 200]);
legend('Simulation','Theory');

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
loglog(f,(Yabs(1:NFFT/2+1)));
ylabel('Y (1/Ohm.m)');xlim([fmin fmax]);
hold on;loglog(fr2,abs(Yor),'r');title('Longitudinal Admittance Y');xlim([fmin fmax]);
legend('Simulation','Theory');
subplot(2,1,2);
semilogx(f,abs(Yangle(1:NFFT/2+1)));
ylabel('\theta Y (Degree)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
hold on;semilogx(fr2,abs(angle(Yor))*(180/pi),'r');xlim([fmin fmax]);ylim([-200 200]);
legend('Simulation','Theory');