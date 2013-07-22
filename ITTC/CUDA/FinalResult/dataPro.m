clear all;
clc;
format long;
load Iters10_32.txt;
y32 = Iters10_32;
load Iters10_48.txt;
y48 = Iters10_48;
load Iters10_128;
y128 = Iters10_128;

x=(0:0.1:1.0);

semilogy(x,y32,':o',x,y48,':x',x,y128,'-.+');

grid on;
title('Iterations = 10');
xlabel('Eb/N0');
ylabel('Ber');
legend({'P=32','P=48','P=128'});
%text(x,y,'picture');
