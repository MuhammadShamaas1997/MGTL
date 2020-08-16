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

for ind=245:256
Bx(ind,:)=0;By(ind,:)=0;Bz(ind,:)=0;
Hx(ind,:)=0;Hy(ind,:)=0;Hz(ind,:)=0;
Dx(ind,:)=0;Dy(ind,:)=0;Dz(ind,:)=0;
Ex(ind,:)=0;Ey(ind,:)=0;Ez(ind,:)=0;
end

Sx=Ex.*conj(Hx);
Sy=Ey.*conj(Hy);
Sz=Ez.*conj(Hz);
S=abs(mag3(Sx,Sy,Sz));

figure;surf(abs(real(Sz(1:200,:))));hold on;
figure;surf(abs(imag(Sz(1:200,:))));hold on;


% hold on;
% plot(0:(2/400):2,(real(Bx(170,:))*B0));
% plot(0:(2/400):2,(imag(Bx(170,:))*B0),'r');
% grid('on');
% xlabel('Distance (mm)');
% ylabel('|Bz (Wb/m^2)|');
% legend('Real','Imaginary')

epsr=1;
a0=1e-1;%0.1mm
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
fmin=1e8;
fmax=1e9;

muinf=1;gamma=.0001;fn=.0001;sigma=-10000;
fr=(1/3000):1e-4:(1/3); 
mur=muinf+(sigma.*fn.*fn)./(-fr.*fr-1i*gamma*fr);
cf=(1+((fr*f0)./(0.2e6)).*((fr*f0)./(0.2e6)));
mur2=mu0+((10000*mu0)./cf)-1i.*cf.*((fr*f0)/(0.2e6));
fr=fr*f0;mur=mur*mu0;sigma=5e-3;eps=(1*eps0)-1i*(sigma./(2*pi*fr));
gammaor=1i.*2.*pi.*fr.*sqrt(mur.*eps);
etaor=sqrt(mur./eps);
Zor=gammaor.*etaor;Yor=gammaor./etaor;
alphaor=real(gammaor);betaor=imag(gammaor);vpor=(2.*pi.*fr)./betaor;


figure;
subplot(2,1,1);loglog(fr,real(mur)/mu0);title('Real \mu');
subplot(2,1,2);loglog(fr,imag(mur)/mu0);title('Imaginary \mu');

figure
plot(abs(By)./abs(Hy))

T=t0/4;
Fs=1/T;
L=256;
obs=90;
L=2^nextpow2(L);
f=Fs/2*linspace(0,1,L/2+1);
%Hy(1:L,obs)=cos(pi*(1:L));
%Hx(1:L,obs)=Hx(1:L,obs)-mean(Hx(1:L,obs));
%Hx(1:L,obs)=exp((-(((1:L)-(L/2)).*((1:L)-(L/2))))/1000000000);
%Hx(1:L,obs)= 1/(4*sqrt(2*pi*0.01))*(exp(-(1:L).^2/(2*0.01)));
%for obs=1:50
FHxi=(fft((Hx(1:L,obs)),L));
FHxo=(fft((Hx(1:L,obs+1)),L));
Gamma=(log(FHxo./FHxi))/(-(1/4)*(a0));
Gamma(1)=0;
%end
% FHxi=mean(FHxis,2);
% FHxo=mean(FHxos,2);
% Gamma=mean(Gammas,2);

% figure;
% subplot(2,1,1)
% hold on;
% plot(t0*(1:L/128),abs(Hx(1:L,obs))*H0);
% xlabel('Time (s)');ylabel('|Hx(t)|')
% subplot(2,1,2)
% hold on;
% plot(t0*(1:L/64),angle(Hx(1:L,obs))*H0);
% xlabel('Time (s)');ylabel('\theta Hx(t)')

figure
subplot(2,1,1)
semilogx(f(1:L/2+1),abs(FHxi(1:L/2+1)))
xlabel('Frequency (s)');ylabel('|Hx (jw)|');xlim([fmin fmax]);
%axis([0 1e9 0 0.2e4])
subplot(2,1,2)
semilogx(f(1:L/2+1),angle(FHxi(1:L/2+1))*(180/pi))
xlabel('Frequency (s)');ylabel('\theta Hx (jw)');xlim([fmin fmax]);
%axis([0 1e9 -200 200])


figure;
subplot(2,1,1)
hold on;
semilogx(f(1:L/2+1),abs(real(Gamma(1:L/2+1))));
%hold on; plot(f(1:L/2+1),2*pi*f(1:L/2+1)/(3e8),'r');
xlabel('Frequency (Hz)');ylabel('\alpha (Np.m^-^1)');xlim([fmin fmax]);
%axis([0 1e9 0 120])

subplot(2,1,2)
semilogx(f(1:L/2+1),abs(imag(Gamma(1:L/2+1))));
%hold on; plot(f(1:L/2+1),2*pi*f(1:L/2+1)/(3e8),'r');
ylabel('\beta (rad.m^-^1)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
%axis([0 1e9 0 120])

figure;hold on;
subplot(2,1,1);semilogx(fr,alphaor);title('alpha');
subplot(2,1,2);semilogx(fr,betaor);title('beta');


figure
f=f';
semilogx(f(1:L/2+1),abs(2*pi*f(1:L/2+1)./(imag(Gamma(1:L/2+1)))));
%hold on; plot(f(1:L/2+1),(3e8),'r');
ylabel('vp (m.s^-^1)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
%axis([0 1e9 0 5e8])

figure;
semilogx(fr,vpor);title('vp');


% FHxi=(fft(Hx(1:L,obs),L));
% FHxo=(fft(Hx(1:L,obs+1),L));
% Gamma=log(FHxo./FHxi)/(-(1/180)*(030*a0));
% Gamma(1)=0;

% X = 1/(4*sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)));

% T=t0;
% Fs=1/T;
% L=256/4;
NFFT=2^nextpow2(L);
f=Fs/2*linspace(0,1,NFFT/2+1);
FEx=fft((Ex(1:L,obs)),NFFT)/L;
FHx=fft((Hx(1:L,obs)),NFFT)/L;
Z=FEx./FHx;
Z=Z*377;
Z(1)=0;

figure;
subplot(2,1,1)
semilogx(f(1:NFFT/2+1),abs(Z(1:NFFT/2+1)));
xlabel('Frequency (Hz)');ylabel('|\eta| (Ohm)');xlim([fmin fmax]);
%axis([1e8 0.5e9 0 100000])

subplot(2,1,2)
semilogx(f(1:NFFT/2+1),(angle(Z(1:NFFT/2+1))*(180/pi)));
ylabel('\Theta \eta (Degree)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
%axis([0 1e9 -200 200])

figure;hold on;
subplot(2,1,1);semilogx(fr,abs(etaor));title('\eta');
subplot(2,1,2);semilogx(fr,angle(etaor)*(180/pi));title('theta \eta');


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
semilogx(f,(Zabs(1:NFFT/2+1)));
ylabel('|Z| (Ohm/m)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
%axis([0 1e9 0 2e5])

subplot(2,1,2);
semilogx(f,(Zangle(1:NFFT/2+1)));
ylabel('\theta Z (Degree)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
%axis([0 1e9 -200 200])

figure;hold on;
subplot(2,1,1);semilogx(fr,abs(Zor));title('|Z|');
subplot(2,1,2);semilogx(fr,angle(Zor)*(180/pi));title('theta Z');

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
semilogx(f,(Yabs(1:NFFT/2+1)));
ylabel('Y (1/Ohm.m)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
%axis([0 1e9 0 0.5])

subplot(2,1,2);
semilogx(f,(Yangle(1:NFFT/2+1)));
ylabel('\theta Y (Degree)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
%axis([0 1e9 -200 200])

figure;hold on;
subplot(2,1,1);semilogx(fr,abs(Yor));title('|Y|');
subplot(2,1,2);semilogx(fr,angle(Yor)*(180/pi));title('theta Y');
