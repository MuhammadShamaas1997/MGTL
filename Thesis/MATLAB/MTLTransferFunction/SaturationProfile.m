% clc;clear all;
% max=2000;
% MMF(1:max+1)=-max/2:max/2;
% K1=1;
% K1exp=0.0002;
% K2=2;
% K2exp=0.0002;
% 
% ind=1;
% 
% for ind=1:max+1
%         Ea(ind)=(150)*(1+exp(-K1exp*MMF(ind))*K1-exp(-K2exp*MMF(ind))*K2);
% end
% %Ea(ind)=150+exp(-0.005*u(1))*150-exp(-0.005*u(1))*300;
% %Ea(ind)=u(1)*0.13-u(1)*u(1)*0.00003
% % plot(1:max+1,Ea)
% 
% 
% for ind=1:max+1
%         muap(ind)=10000*exp(-(abs(MMF(ind)-250))*(abs(MMF(ind)-250))/50000);
%         muan(ind)=-10000*exp(-(abs(MMF(ind)+250))*(abs(MMF(ind)+250))/50000);
% end
% 
% hold on;
% plot((-max/2):1:(max/2),muap)
% %plot((-max/2):1:(max/2),muan,'r')
% grid on;
% xlabel('H (A/m)');
% ylabel('Magnetic Susceptibility \chi_m');
% %legend('\chi_m^+','\chi_m^-')

% B(1)=0;
% Hmax=1000;
% Happ=[0:Hmax/4 Hmax:-1:(-Hmax) (-Hmax):1:(Hmax)];
% hold on;
% for ind=2:length(Happ)
%     if (Happ(ind)>Happ(ind-1))
%         munet=(muap(Happ(ind)+(Hmax+1))+1);
%         dB=(Happ(ind)-Happ(ind-1))*munet*(4*pi*1e-7);
%     else
%         munet=(muan(Happ(ind)+(Hmax+1))+1);
%         dB=(-Happ(ind)+Happ(ind-1))*munet*(4*pi*1e-7);
%     end
%     B(ind)=B(ind-1)+dB;
% end
% 
% plot(Happ,B);
% grid('on');
% xlabel('H (A/m)');
% ylabel('B (T)');

% subplot(2,3,1)
% plot(ScopeData5.time,ScopeData5.signals(1,1).values)
% xlabel('Time (s)');ylabel('Speed (rpm)');%axis([0 0.04 0 2500]);
% subplot(2,3,2)
% plot(ScopeData5.time,ScopeData5.signals(1,2).values)
% xlabel('Time (s)');ylabel('IFshunt (A)');%axis([0 0.04 2 5]);
% subplot(2,3,3)
% plot(ScopeData5.time,ScopeData5.signals(1,3).values)
% xlabel('Time (s)');ylabel('IFSeries (A)');%axis([0 0.04 20 65]);
% subplot(2,3,4)
% plot(ScopeData5.time,ScopeData5.signals(1,4).values)
% xlabel('Time (s)');ylabel('Net MMF(A.t)');%axis([0 0.04 2000 6500]);
% subplot(2,3,5)
% plot(ScopeData5.time,ScopeData5.signals(1,5).values)
% xlabel('Time (s)');ylabel('Net EMF(V)');%axis([0 0.04 50 175]);
% subplot(2,3,6)
% plot(ScopeData5.time,ScopeData5.signals(1,6).values)
% xlabel('Time (s)');ylabel('Terminal Voltage (V)');%axis([0 0.04 50 150]);

subplot(2,1,1)
mmf=ScopeData5.signals(1,1).values;
emf=ScopeData5.signals(1,5).values;
plot(mmf,emf);
ylabel('Generated EMF (V)');xlabel('Net MMF (A)');
%axis([400 2200 0 150]);
subplot(2,1,2)
plot(ScopeData5.time,ScopeData5.signals(1,5).values)
xlabel('Time (s)');ylabel('Terminal Voltage (V)');%axis([0 0.04 50 150]);
