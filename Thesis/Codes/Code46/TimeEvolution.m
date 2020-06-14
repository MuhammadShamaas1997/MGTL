clc;clear all;
f=fopen('TimeEvolution.txt');
l=fgetl(f);
i=1;
while ischar(l)
    %%disp(l);
    text{i}=l;
    data{i}=sscanf(text{i},'%f , %f , %f , %f , %f , %f , %f , %f');
    A1(i)=data{i}(1);
    A2(i)=data{i}(2);
    A3(i)=data{i}(3);
    A4(i)=data{i}(4);
    A5(i)=data{i}(5);
    A6(i)=data{i}(6);
    A7(i)=data{i}(7);
    A8(i)=data{i}(8);
    l=fgetl(f);
    i=i+1;
end

hold on;

subplot(2,2,1)
hold on;plot(mag(A1,A2));%Im
ylabel('Im (V)');
xlabel('Time')

subplot(2,2,2)
hold on;plot(mag(A3,A4));%Vm
ylabel('Vm (A)');
xlabel('Time')

subplot(2,2,3)
hold on;plot(mag(A5,A6));%Ie
ylabel('Ie (A)');
xlabel('Time')

subplot(2,2,4)
hold on;plot(mag(A7,A8));%Ve
ylabel('Ve (V)');
xlabel('Time')

