subplot(4,1,1)
semilogx(f,2*abs(YHz(1:NFFT/2+1)));
% plot(f,2*angle(YHz(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Hz(f)|');
%axis([])

subplot(4,1,2)
semilogx(f,2*abs(YEz(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Ez(f)|');

subplot(4,1,3)
semilogx(f,2*abs(Z(1:NFFT/2+1)));
xlabel('Frequency (Hz)')
ylabel('|Z(f)|');
%axis([0 2e9 0  10e6])

subplot(4,1,4)
semilogx(f,2*angle(Z(1:NFFT/2+1))*(180/pi),'.-');
ylabel('\Theta Z(f)');
xlabel('Frequency (Hz)')
%axis([0 2e9 -200  200])