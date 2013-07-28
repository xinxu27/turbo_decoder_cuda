clear all;
clc;
format long;
%load Iters10_1;
%y1 = Iters10_1;
load Iter32Eb01.txt;
y32 = Iter32Eb01;
load Iter48Eb01.txt;
y48 = Iter48Eb01;
load Iter64Eb01.txt;
y64 = Iter64Eb01;
load Iter96Eb01.txt;
y96 = Iter96Eb01;
load Iter128Eb01.txt;
y128 = Iter128Eb01;


x=(1:1:15);

semilogy(x,y32,'-ok',x,y48,':dr',x,y64,'--p',x,y96,'-vm',x,y128,'-.+','linewidth',2);

grid on;
title('Eb/N0 = 0.1db');
xlabel('Iteration');
ylabel('BER');
legend({'P=32','P=48','P=64','P=96','P=128'});
%text(x,y,'picture');
