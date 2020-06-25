clc;clear all;
%hold on;
%%
f=fopen('F1.txt');
File=1;
l=fgetl(f);
in=1;
Prev=1600*(File-1);
while ischar(l)
%for kj=1:81
    %%disp(l);
    text{in}=l;
    data{in}=sscanf(text{in},'%f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f');
    Hx(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(3)+1i*data{in}(4);
    Hy(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(5)+1i*data{in}(6);
    Hz(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(7)+1i*data{in}(8);
    
    Bx(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(9)+1i*data{in}(10);
    By(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(11)+1i*data{in}(12);
    Bz(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(13)+1i*data{in}(14);
    
    Ex(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(15)+1i*data{in}(16);
    Ey(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(17)+1i*data{in}(18);
    Ez(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(19)+1i*data{in}(20);
    
    Dx(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(21)+1i*data{in}(22);
    Dy(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(23)+1i*data{in}(24);
    Dz(data{in}(1)-Prev,(((data{in}(2)+10)*4)+1))=data{in}(25)+1i*data{in}(26);
    
    l=fgetl(f);
    in=in+1;
end


epsr=10;
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
sigmaD0=(epsr*eps0*c0)/a0;

%subplot(9,1,1);
% hold on;
% for tr=1:5:100
% plot(0:0.025:2,log((abs(Bz(tr,:))+1e-59)*B0))
% pause(0.01);
% end
% xlabel('Length (mm)')
% ylabel('log |Bz|')
% grid('on')

% hold on;
% plot((1:1600)*t0,log((abs(Bz(1:1600,40)))*B0))
% xlabel('Time (s)')
% ylabel('log |Bz|')
% grid('on')



%hold on;

% subplot(2,1,1)
% hold on;%plot(mag3(mag(A1,A2),mag(A3,A4),mag(A5,A6)));%H
% for m=1:length(A3)
%         Hy(m)=A3(m)+i*A4(m);
%         %Hz=(mag3(mag(A1,A2),mag(A3,A4),mag(A5,A6)));
%         Hy(m)=Hy(m)*H0;
% end


T=t0;
Fs=1/T;
L=45;
NFFT=2^nextpow2(L);
FEx=fft(Ex(1:L,1)*E0,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);




%plot(f,2*abs(FEx(1:NFFT/2+1)));
% plot(f,2*angle(FEx(1:NFFT/2+1)));
% xlabel('Frequency (Hz)')
% ylabel('|Ex(f)|');
%axis([2e9 15e9 0  1e-3])

%subplot(2,2,2)
%hold on;%plot(mag3(mag(A7,A8),mag(A9,A10),mag(A11,A12)));%B

% subplot(2,1,2)
% hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
% for m=1:length(A15)
%         Ey(m)=A15(m)+i*A16(m);
%         Ey(m)=Ey(m)*E0;
% end


FHx=fft(Hx(1:L,1)*H0,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);

% plot(f,2*abs(FHx(1:NFFT/2+1)));
% xlabel('Frequency (Hz)')
% ylabel('|Hx(f)|');
%axis([2e9 15e9 0  1e3])

Z=FEx./FHx;

subplot(2,1,1)
hold on;%plot(mag3(mag(A13,A14),mag(A15,A16),mag(A17,A18)));%E
plot(f,2*abs(Z(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Z_w (Ohm)|');
%axis([0 10e9 5.5e4  5.7e4])
axis([0 1e11 0  2e5])

% subplot(2,2,4)
% hold on;plot(mag3(mag(A19,A20),mag(A21,A22),mag(A23,A24)));%D

subplot(2,1,2)
plot(f,angle(Z(1:NFFT/2+1))*(180/pi),'.-');
ylabel('\Theta Z_w (f)');
xlabel('Frequency (Hz)')
%axis([0 10e9 -200 0])
axis([0 1e11 -200 200])
% 
% % %smithchart(z2gamma(Z,120*pi*22))
% 
% 
Lm=(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1))./(2*pi*(f+1)');
G=real(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1));
Cm=imag(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1))./(2*pi*(f+1)');

subplot(3,1,1)
plot(f,abs(Lm(1:NFFT/2+1)));
ylabel('Lm(H/m)');
xlabel('Frequency (Hz)')
axis([0 1e11 0 2e9])

subplot(3,1,2)
plot(f,abs(G(1:NFFT/2+1)));
ylabel('G(S/m)');
xlabel('Frequency (Hz)')
axis([0 1e11 0 15])

subplot(3,1,3)
plot(f,abs(Cm(1:NFFT/2+1)));
ylabel('Cm(F/m)');
xlabel('Frequency (Hz)')
axis([0 1e11 0 0.05])