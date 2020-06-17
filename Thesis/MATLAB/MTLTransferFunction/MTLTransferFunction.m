L=0.5;
G=0.4;
C=0.1;
H=tf([L*C L*G 0],[1]);
%bode(G)
%nyquist(G)


Y=tf([1],[L*C L*G 0]);
bode(Y)