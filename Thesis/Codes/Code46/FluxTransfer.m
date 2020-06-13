clc;clear all;
f=fopen('Flux.txt');
l=fgetl(f);
i=1;
while ischar(l)
    %%disp(l);
    text{i}=l;
    data{i}=sscanf(text{i},'%f , %f , %f');
    A1(i)=data{i}(1);
    A2(i)=data{i}(2);
    A3(i)=data{i}(3);
    l=fgetl(f);
    i=i+1;
end

hold on;

subplot(3,1,1)
A=abs(A2)/max(abs(A2));%Flux_in
B=abs(A3)/max(abs(A2));%Flux_out

semilogy(A1,-A2);
axis([2.5e-3 2.5e-2 0.01 10]);
ylabel('log |Sin|');
xlabel('frequency');


subplot(3,1,2)
semilogy(A1,-A3);
axis([2.5e-3 2.5e-2 0.0000001 0.01]);
ylabel('log |Sout|');
xlabel('frequency');


subplot(3,1,3)
C=B./A;
for i=1:length(C)
    if (C(i)>1)
        %C(i)=0;
    end
end
semilogy(A1,C);
axis([2.5e-3 2.5e-2 0.00001 0.001]);
ylabel('log |Sout/Sin|');
xlabel('Length (mm)');
