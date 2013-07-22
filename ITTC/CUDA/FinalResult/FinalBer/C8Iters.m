clear all;
clc;
format long;
load Iters8_32;
y32 = Iters8_32;
load Iters8_48;
y48 = Iters8_48;
load Iters8_64;
y64 = Iters8_64;
load Iters8_96;
y96 = Iters8_96;
load Iters8_128;
y128 = Iters8_128;


x=(0:0.1:1.0);

semilogy(x,y32,'-ok',x,y48,':x',x,y64,'--p',x,y96,'-*',x,y128,'-.+','linewidth',2);

grid on;
title('Iterations = 8');
xlabel('Eb/N0');
ylabel('Ber');
legend({'P=32','P=48','P=64','P=96','P=128'});
%text(x,y,'picture');
