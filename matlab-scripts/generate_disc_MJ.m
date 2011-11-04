clear

filename = 'disc_200nm_M.init';
filenameJ = 'disc_200nm_J.init';
Ms = 8.6e5;
dia = 200e-9;
diaJ = 50e-9;
I = 1;
dx = 5e-9;

r = dia/2;
Nr = round(r/dx);
Nx = 2*Nr+3;
Ny = Nx;
x = [(-Nx+1)/2:(Nx-1)/2];
y = [(-Ny+1)/2:(Ny-1)/2];
[X Y] = meshgrid(x,y);

A = zeros(Ny,Nx);
A = (X.^2+Y.^2 <= Nr^2);

M = zeros(Ny,Nx,3);

for r = [1:Ny]
    for c = [1:Nx]
        if A(r,c) == 1
            theta = pi/2;
            phi = pi/4;
            theta = pi*rand;
            phi = 2*pi*rand;
            M(r,c,:) = Ms * [sin(theta) * cos(phi)
                             sin(theta) * sin(phi)
                             cos(theta)];
        end
    end
end

J = zeros(Ny,Nx,3);
rJ = diaJ/2;
NrJ = round(rJ/dx);
currentDensity = I/(pi*rJ^2);
J = currentDensity * (X.^2+Y.^2 <= NrJ^2);
Jfile = fopen(filenameJ, 'w');
fprintf(Jfile, 'Nx=%d\n', Nx);
fprintf(Jfile, 'Ny=%d\n', Ny);
fprintf(Jfile, '\n# start field\n', dx);
for r = [1:Ny]
    for c = [1:Nx]
        fprintf(Jfile, '%g ', J(r,c));
    end
    fprintf(Jfile, '\n');
end
fprintf(Jfile, '# end field\n');
fclose(Jfile);
disp([filenameJ, ' file written'])


Mfile = fopen(filename, 'w');
fprintf(Mfile, 'Nx=%d\n', Nx);
fprintf(Mfile, 'Ny=%d\n', Ny);
% fprintf(Mfile, 'cellSize=%g\n', dx);
fprintf(Mfile, '\n# start field\n', dx);

for r = [1:Ny]
    for c = [1:Nx]
        m = M(r,c,:);
        fprintf(Mfile, '%g %g %g\n', m(1), m(2), m(3));
    end
    fprintf(Mfile, '\n');
end

fprintf(Mfile, '# end field\n');
fclose(Mfile);
disp([filename, ' file written'])



imagesc(J+A*max(J(:))); axis image;

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
