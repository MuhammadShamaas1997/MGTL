w=1e-5;
t0=5;
t=0:0.001:10;
g=exp(-1i*w*t-((t-t0).*(t-t0))./(2*w*w));
hold on;
plot(t,abs(g));
plot(t,real(g),'r')
plot(t,imag(g),'c')