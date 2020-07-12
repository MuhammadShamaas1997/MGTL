clc;clear all;
Nf=1000;
If=0:0.5:7;
Ea=[0 20 40 60 80 95 107 115 121 125 129 133 136 139.5 142];
plot(Nf*If,Ea);
xlabel('Field MMF (A.t)');ylabel('Generated EMF (V)');title('Magnetization Curve of Generator at 1800rpm');
K=1;
Omega=1800;
phi=Ea/(Omega*K);
% subplot(2,1,2)
% plot(Nf*If,phi);
% xlabel('Field MMF (A.t)');ylabel('Flux (V.s)');title('Flux of Generator at 1800rpm');
for ind=1:14
Rm(ind)=(Nf*(If(ind+1)-If(ind)))./(phi(ind+1)-phi(ind));
end
figure;
plot(Nf*If(1:14),Rm)
hold on;
xlabel('Field MMF (A.t)');ylabel('Magnetic Reluctance (H^-^1)');title('Reluctance and cubic fitting function at 1800rpm');
% -1.9e6x*x*x+0.024*x*x-30x+4.8e4