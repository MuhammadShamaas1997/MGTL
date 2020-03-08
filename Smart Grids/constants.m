clc;clear;
Grid_Angle_init = pi;

Pb=2.3e6;
Vp=(690*(2^0.5))/(3^0.5);
Fb=60;

Vbase=Vp/(2^.5);
Ibase=(Pb/3)/(690/(3^0.5));

Lgrid=0.1098e-3;

%switching frequency
Fs=1e3;

% sampling frqquency for controller
Fsampling=2e3;

Kp=1;
Ki=2;