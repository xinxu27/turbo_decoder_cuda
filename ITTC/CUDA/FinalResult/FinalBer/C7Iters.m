clear all;
clc;
format long;
%load Iters7_1;
%y1 = Iters7_1;
load Iters7_32;
y32 = Iters7_32;
load Iters7_48;
y48 = Iters7_48;
load Iters7_64;
y64 = Iters7_64;
load Iters7_96;
y96 = Iters7_96;
load Iters7_128;
y128 = Iters7_128;


x=(0:0.1:1.0);

semilogy(x,y32,'-ok',x,y48,':dr',x,y64,'--p',x,y96,'-vm',x,y128,'-.+','linewidth',2);

grid on;
title('Iteration = 7');
xlabel('Eb/N0');
ylabel('BER');
legend({'P=32','P=48','P=64','P=96','P=128'});
%text(x,y,'picture');
