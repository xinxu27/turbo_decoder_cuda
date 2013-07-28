clear all;
clc;
format long;
%load Iters6_1;
%y1 = Iters6_1;
load Iters6_32;
y32 = Iters6_32;
load Iters6_48;
y48 = Iters6_48;
load Iters6_64;
y64 = Iters6_64;
load Iters6_96;
y96 = Iters6_96;
load Iters6_128;
y128 = Iters6_128;


x=(0:0.1:1.0);

semilogy(x,y32,'-ok',x,y48,':dr',x,y64,'--p',x,y96,'-vm',x,y128,'-.+','linewidth',2);

grid on;
title('Iteration = 6');
xlabel('Eb/N0');
ylabel('BER');
legend({'P=32','P=48','P=64','P=96','P=128'});
%text(x,y,'picture');
