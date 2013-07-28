clear all;
clc;
format long;
%load Iters10_1;
%y1 = Iters10_1;
load Iter32Eb06.txt;
y32 = Iter32Eb06;
load Iter48Eb06.txt;
y48 = Iter48Eb06;
load Iter64Eb06.txt;
y64 = Iter64Eb06;
load Iter96Eb06.txt;
y96 = Iter96Eb06;
load Iter128Eb06.txt;
y128 = Iter128Eb06;


x=(1:1:15);

semilogy(x,y32,'-ok',x,y48,':dr',x,y64,'--p',x,y96,'-vm',x,y128,'-.+','linewidth',2);

grid on;
title('Eb/N0 = 0.6db');
xlabel('Iteration');
ylabel('BER');
legend({'P=32','P=48','P=64','P=96','P=128'});
%text(x,y,'picture');
