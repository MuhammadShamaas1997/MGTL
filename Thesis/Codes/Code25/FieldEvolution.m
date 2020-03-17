clc;clear all;
f=fopen('FieldEvolutionIn.txt');
l=fgetl(f);
i=1;
while ischar(l)
    %%disp(l);
    text{i}=l;
    data{i}=sscanf(text{i},'%f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f');
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
    A13(i)=data{i}(13);
    A14(i)=data{i}(14);
    A15(i)=data{i}(15);
    A16(i)=data{i}(16);
    A17(i)=data{i}(17);
    A18(i)=data{i}(18);
    A19(i)=data{i}(19);
    A20(i)=data{i}(20);
    A21(i)=data{i}(21);
    A22(i)=data{i}(22);
    A23(i)=data{i}(23);
    A24(i)=data{i}(24);
    l=fgetl(f);
    i=i+1;
end

hold on;

subplot(4,3,1)
hold on;plot(mag(A1,A2));%Hx
subplot(4,3,2)
hold on;plot(mag(A3,A4));%Hy
subplot(4,3,3)
hold on;plot(mag(A5,A6));%Hz
subplot(4,3,4)
hold on;
% plot(mag(A7,A8)./mag(A1,A2));%Bx
plot((A7)./(A1),'o');
plot((A8)./(A2),'*');
subplot(4,3,5)
hold on;plot(mag(A9,A10)./mag(A3,A4));%By
subplot(4,3,6)
hold on;plot(mag(A11,A12)./mag(A5,A6));%Bz
subplot(4,3,7)
hold on;plot(mag(A13,A14));%Ex
subplot(4,3,8)
hold on;plot(mag(A15,A16));%Ey
subplot(4,3,9)
hold on;plot(mag(A17,A18));%Ez
subplot(4,3,10)
hold on;plot(mag(A19,A20)./mag(A13,A14));%Dx
subplot(4,3,11)
hold on;plot(mag(A21,A22)./mag(A15,A16));%Dy
subplot(4,3,12)
hold on;plot(mag(A23,A24)./mag(A17,A18));%Dz
