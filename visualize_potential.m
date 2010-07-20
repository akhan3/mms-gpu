clear
load charge.dat
load potential.dat

[ylen, xlen] = size(charge);
x = 0:xlen-1;
y = 0:ylen-1;
[X, Y] = meshgrid(x,y);

subplot(121)
imagesc(x,y,charge); axis image xy;
subplot(122)
imagesc(x,y,potential); axis image xy;
