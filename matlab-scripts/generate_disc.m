clear

A = [ 
 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 1 1 1 1 1 0 0 0 0
 0 0 0 1 1 1 1 1 1 1 0 0 0
 0 0 1 1 1 1 1 1 1 1 1 0 0
 0 1 1 1 1 1 1 1 1 1 1 1 0
 0 1 1 1 1 1 1 1 1 1 1 1 0
 0 1 1 1 1 1 1 1 1 1 1 1 0
 0 1 1 1 1 1 1 1 1 1 1 1 0
 0 1 1 1 1 1 1 1 1 1 1 1 0
 0 0 1 1 1 1 1 1 1 1 1 0 0
 0 0 0 1 1 1 1 1 1 1 0 0 0
 0 0 0 0 1 1 1 1 1 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0
];

R = 13;
C = 13;

M = zeros(13,13,3);
r = 5;
x = [(-length(M)+1)/2:(length(M)-1)/2];
y = x;
[X Y] = meshgrid(x,y);

Ms = 8.6e5;

for r = [1:R]
    for c = [1:C]
        if A(r,c) == 1
            theta = pi/2;
            phi = pi/3;
%             theta = pi*rand;
%             phi = 2*pi*rand;
            M(r,c,:) = Ms * [sin(theta) * cos(phi)
                             sin(theta) * sin(phi)
                             cos(theta)];
        end
    end
end

filename = 'disc_13x13.init';
Mfile = fopen(filename, 'w');
for r = [1:R]
    for c = [1:C]
        m = M(r,c,:);
        fprintf(Mfile, '%g %g %g\n', m(1), m(2), m(3));
    end
    fprintf(Mfile, '\n');
end
disp([filename, ' file written'])

return

M = 2*double(sqrt(X.^2+Y.^2) == r);
M = M+double(sqrt(X.^2+Y.^2) < r);
% [C h] = contourf(x,y,M);
% text_handle = clabel(C,h);

% figure
imagesc(M)
axis image 
axis xy
grid on
colormap('spring')
%colorbar
