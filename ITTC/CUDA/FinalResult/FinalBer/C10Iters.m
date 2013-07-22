clear all;
clc;
format long;
load Iters10_32;
y32 = Iters10_32;
load Iters10_48;
y48 = Iters10_48;
load Iters10_64;
y64 = Iters10_64;
load Iters10_96;
y96 = Iters10_96;
load Iters10_128;
y128 = Iters10_128;


x=(0:0.1:1.0);

semilogy(x,y32,'-ok',x,y48,':x',x,y64,'--p',x,y96,'-*',x,y128,'-.+','linewidth',2);

grid on;
title('Iterations = 8');
xlabel('Eb/N0');
ylabel('Ber');
legend({'P=32','P=48','P=64','P=96','P=128'});
%text(x,y,'picture');
