clc
clear
F0=1;
Grid_Angle_init = pi;
iq_ref_pu = 0;
id_ref_pu = sqrt(2);
% rectifer or inverter selection
% "1" --- rectifier
% "0" --- inverter
Rectifier_or_inverter=0;


Pb=2.3e6;% active power at the PCC
Vp=(690*(2^0.5))/(3^0.5);%563.3826
Vbase=Vp/(2^.5);
Fb=60;
Ibase=(Pb/3)/(690/(3^0.5));
Zb=Vbase/Ibase;
Omega_b=2*pi*Fb;
XLb=Zb;
XCb=Zb;

% grid-side filter inductor
Lgrid=0.1098e-3;
XLgrid=0.2*XLb;
XLgrid_pu=XLgrid/Zb;
Rgrid_pu=0.1;
Rgrid=Rgrid_pu*Zb;
Power_loss=Rgrid_pu*Pb;

% %dc link
% XCdc=0.3*XCb;
% Cdc=1/(Omega_b*XCdc);
% Vdc_ref=3.062;% in per unit
% Vdc_capacitor_initial=Vdc_ref*Vbase;
% 
% %dc supply voltage calculation according to series resistor
% Rdc_series=0.1*Zb; % might be changed here
% 
% if Rectifier_or_inverter==1 % for rectifier operation case
%     Idoc=(Pb - Power_loss)/(Vdc_ref*Vbase);
%     Vdoc_supply=Vdc_ref*Vbase - Idoc*Rdc_series;
% elseif Rectifier_or_inverter==0% for inverter operation case
%     Idoc=(Pb + Power_loss)/(Vdc_ref*Vbase);
%     Vdoc_supply=Vdc_ref*Vbase + Idoc*Rdc_series;
% end

%switching frequency
Fs=2040;% 34 times fundamental frequency

% sampling frqquency for controller
Fsampling=10e3;% at least 2 times Fs recommended


% Kp=1;
% Ki=2;
Kp=1;
Ki=2;