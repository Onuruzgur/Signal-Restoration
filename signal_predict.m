clc;  clear all;  close all;

data     = load('SigRes.mat');
H        = data.H;
y        = data.y;
x        = data.x;
L        = data.L;
varyans=25*10^-4;
Xml=inv((H'\H))*(H'*y);
figure();
plot(x);
title("X observation");
figure();
plot(y);
title("y observation");
t=linspace(1,256,256);
figure();
plot(Xml,'-b')
hold on
plot(x,'-.r')
title("Xml and truth x");
lambda=[1,10,100,1000];

Lt=L';
b=(H'*y);
Ex1=inv(lambda(1)*Lt*L);
Ex2=inv(lambda(2)*Lt*L);
Ex3=inv(lambda(3)*Lt*L);
Ex4=inv(lambda(4)*Lt*L);
Ex1i=(lambda(1)*Lt.*L);
Ex2i=(lambda(2)*Lt.*L);
Ex3i=(lambda(3)*Lt.*L);
Ex4i=(lambda(4)*Lt.*L);

Xmap1=inv(H'*H+ varyans*Ex1i)*(b);
Xmap2=inv(H'*H+ varyans*Ex2i)*(b);
Xmap3=inv(H'*H+ varyans*Ex3i)*(b);
Xmap4=inv(H'*H+ varyans*Ex4i)*(b);

figure()
plot(Xmap1,'-b')
hold on
plot(x,'-.r')
title("Xmap with 𝜆=1");
figure()
plot(Xmap2,'-b')
hold on
plot(x,'-.r')
title("Xmap with 𝜆=10");
figure()
plot(Xmap3,'-b')
hold on
plot(x,'-.r')
title("Xmap with 𝜆=100");
figure()
plot(Xmap4,'-b')
hold on
plot(x,'-.r')
title("Xmap with 𝜆=1000");
mean=0;
I=eye(256);
temp=inv(H*Ex1*H'+varyans*I);
Xmmse1=Ex1*H'*temp*y;
Xmmse2=Ex2*H'*temp*y;
Xmmse3=Ex3*H'*temp*y;
Xmmse4=Ex4*H'*temp*y;

figure()
plot(Xmmse1,'-b')
hold on
plot(x,'-.r')
title("Xmmse with 𝜆=1");
figure()
plot(Xmmse2,'-b')
hold on
plot(x,'-.r')
title("Xmmse with 𝜆=10");
figure()
plot(Xmmse3,'-b')
hold on
plot(x,'-.r')
title("Xmmse with 𝜆=100");
figure()
plot(Xmmse4,'-b')
hold on
plot(x,'-.r')
title("Xmmse with 𝜆=1000");


PSNR_Xml=20*log10((max(x)*sqrt(256))/(norm(x-Xml)));

PSNR_Xmap1=20*log10((max(x)*sqrt(256))/(norm(x-Xmap1)));
PSNR_Xmap2=20*log10((max(x)*sqrt(256))/(norm(x-Xmap2)));
PSNR_Xmap3=20*log10((max(x)*sqrt(256))/(norm(x-Xmap3)));
PSNR_Xmap4=20*log10((max(x)*sqrt(256))/(norm(x-Xmap4)));

PSNR_Xmmse1=20*log10((max(x)*sqrt(256))/(norm(x-Xmmse1)));
PSNR_Xmmse2=20*log10((max(x)*sqrt(256))/(norm(x-Xmmse2)));
PSNR_Xmmse3=20*log10((max(x)*sqrt(256))/(norm(x-Xmmse3)));
PSNR_Xmmse4=20*log10((max(x)*sqrt(256))/(norm(x-Xmmse4)));


fprintf('PSNR_Xml=%.4f\n',PSNR_Xml);
fprintf('PSNR_Xmap(𝜆=1)=%.4f\n',PSNR_Xmap1);
fprintf('PSNR_Xmap(𝜆=10)=%.4f\n',PSNR_Xmap2);
fprintf('PSNR_Xmap(𝜆=100)=%.4f\n',PSNR_Xmap3);
fprintf('PSNR_Xmap(𝜆=1000)=%.4f\n',PSNR_Xmap4);

fprintf('PSNR_Xmmse(𝜆=1)=%.4f\n',PSNR_Xmmse1);
fprintf('PSNR_Xmmse(𝜆=10)=%.4f\n',PSNR_Xmmse2);
fprintf('PSNR_Xmmse(𝜆=100)=%.4f\n',PSNR_Xmmse3);
fprintf('PSNR_Xmmse(𝜆=1000)=%.4f\n',PSNR_Xmmse4);