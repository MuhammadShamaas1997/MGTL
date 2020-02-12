clc;clear all;
f=fopen('FieldEvolution.txt');
l=fgetl(f);
i=1;
while ischar(l)
    %%disp(l);
    text{i}=l;
    data{i}=sscanf(text{i},'%f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f');
    A1(i)=data{i}(1);
    A2(i)=data{i}(2);
    A3(i)=data{i}(3);
    A4(i)=data{i}(4);
    A5(i)=data{i}(5);
    A6(i)=data{i}(6);
    A7(i)=data{i}(7);
    A8(i)=data{i}(8);
    A9(i)=data{i}(9);
    A10(i)=data{i}(10);
    A11(i)=data{i}(11);
    A12(i)=data{i}(12);
    l=fgetl(f);
    i=i+1;
end

hold on;

subplot(2,3,1)
hold on;plot(mag(A1,A2));%Hx
subplot(2,3,2)
hold on;plot(mag(A3,A4));%Hy
subplot(2,3,3)
hold on;plot(mag(A5,A6));%Hz
subplot(2,3,4)
hold on;plot(mag(A7,A8));%Bx
subplot(2,3,5)
hold on;plot(mag(A9,A10));%By
subplot(2,3,6)
hold on;plot(mag(A11,A12));%Bz
