clear all;
clc;
format long;
%load Iters11_1;
%y1 = Iters11_1;
load Iters11_32;
y32 = Iters11_32;
load Iters11_48;
y48 = Iters11_48;
load Iters11_64;
y64 = Iters11_64;
load Iters11_96;
y96 = Iters11_96;
load Iters11_128;
y128 = Iters11_128;


x=(0:0.1:1.0);

semilogy(x,y32,'-ok',x,y48,':dr',x,y64,'--p',x,y96,'-vm',x,y128,'-.+','linewidth',2);

grid on;
title('Iteration = 10');
xlabel('Eb/N0');
ylabel('BER');
legend({'P=32','P=48','P=64','P=96','P=128'});
%text(x,y,'picture');
