clear all;
clc;
format long;
%load Iters12_1;
%y1 = Iters12_1;
load Iters12_32;
y32 = Iters12_32;
load Iters12_48;
y48 = Iters12_48;
load Iters12_64;
y64 = Iters12_64;
load Iters12_96;
y96 = Iters12_96;
load Iters12_128;
y128 = Iters12_128;


x=(0:0.1:1.0);

semilogy(x,y32,'-ok',x,y48,':dr',x,y64,'--p',x,y96,'-vm',x,y128,'-.+','linewidth',2);

grid on;
title('Iteration = 12');
xlabel('Eb/N0');
ylabel('BER');
legend({'P=32','P=48','P=64','P=96','P=128'});
%text(x,y,'picture');
