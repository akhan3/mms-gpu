clear

filename = 'disc_500nm_M_sd.init';
filenameJ = 'disc_500nm_J_50nm.init';
Ms = 8.6e5;
dia = 500e-9;
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
            random = 0;
            if(random)
                theta = pi*rand;
                phi = 2*pi*rand;
            else
                theta = pi/2;
                phi = pi/4;
            end


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


%% subsample
sf = 3;
x = x(1:sf:end);
y = y(1:sf:end);
[X,Y] = meshgrid(x,y);
Z = 0*X;
Nx = length(x);
Ny = length(y);
Mx = M(1:sf:end, 1:sf:end, 1);
My = M(1:sf:end, 1:sf:end, 2);
Mz = M(1:sf:end, 1:sf:end, 3);

%% plot
clf
q1 = quiver3(X,Y,Z, Mx(:,:,1), My(:,:,1), Mz(:,:,1), .5);
%     q1 = quiver3(X,Y,Z, Mx(:,:,1), My(:,:,1), zeros(size(Mz(:,:,1))), .5);
set(q1,'color','black','linewidth',2);
grid off;
xlabel('x'); ylabel('y'); title('Magnetization (\phi - angle in plane)');
view(3);
set(gca,'visible','off')
hold on
phi = atan2(My,Mx) * 180/pi;
q3 = pcolor(X(1,:),Y(:,1), phi(:,:,1) + J(1:sf:end,1:sf:end));
hold off;
shading interp;
axis equal tight xy;
grid off;
colormap(hsv); 
% colorbar
caxis([-180,180])    
hold off

return
