clc;clear all;
f=fopen('AG.MAT');
l=fgetl(f);
l=fgetl(f);
l=fgetl(f);
d=sscanf(l,'%s%d');
num=d(7);
l=fgetl(f);
i=1;
for i=1:num
    %disp(l);
    text{i}=l;
    data{i}=sscanf(text{i},' %s %f %f %f %f');
    w1(i)=data{i}(7);
    n(i)=data{i}(8);
    k(i)=data{i}(9);
    eps(i) = ((n(i)+1i*k(i))*(n(i)+1i*k(i)));
    freq(i) = (2.998e8 / (w1(i)*1e-9));
    l=fgetl(f);
    i=i+1;
end

% wl, n, k = np.array(data).T
% eps = ((n+1j*k)**2)[::-1]
% freq = (2.998e8 / (wl*1e-9))[::-1]

hold on;
plot(freq,imag(eps))
plot(freq,real(eps),'r')
plot(freq,abs(eps),'k')

