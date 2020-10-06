hold on;
L2=32768*4;
L=L2;
L=2^nextpow2(L);
f=Fs/2*linspace(0,1,L/2+1);
plot(f(1:L/2+1),abs(imag(Gamma(1:L/2+1))));
% ylabel('\beta (rad.m^-^1)');xlim([fmin fmax]);

