L=0.5;
G=0.4;
C=0.1;

num={1 [-C -G] ; 0 1};
den={1};
X=tf(num,den);
%stepplot(Y)

num={1 0 ; [-L 0] 1};
den={1};
Y=tf(num,den);

sys=[X;Y;X];

num={1 0 ; [-C -G] 1};
den={1};
Y=tf(num,den);
bode(Y)
