clear

load dynamics.dat;
    load dynamics.dat
    tindex = dynamics(:,1);
    time = dynamics(:,2);
    dt = dynamics(:,3);
    energy = dynamics(:,4);
    dE = dynamics(:,5); dE(1) = dE(2);
    Mx = dynamics(:,6);
    My = dynamics(:,7);
    Mz = dynamics(:,8);
    M  = dynamics(:,9);
    torque_max  = dynamics(:,10);

    tdim = length(time);
    zdim = 1;
    ydim = 16;
    xdim = ydim;
    x = 0:xdim-1;
    y = 0:ydim-1;
    z = 0:zdim-1;
    [X,Y,Z] = meshgrid(x,y,z);
    required_lines = ydim*xdim*zdim*tdim


M_yxzt = load('Mdynamics.dat');
M_yxzt = reshape(M_yxzt', 3,ydim,xdim,zdim, tdim);
Mx = shiftdim(M_yxzt(1,:,:,:,:), 1);
My = shiftdim(M_yxzt(2,:,:,:,:), 1);
Mz = shiftdim(M_yxzt(3,:,:,:,:), 1);


figure;
set(gcf, 'OuterPosition', [0 0 1280 800]);
subplot(121);
    q1 = quiver3(X,Y,Z, Mx(:,:,:,1), My(:,:,:,1), Mz(:,:,:,1), .5);
    axis tight equal;
    xlabel('x'); ylabel('y'); zlabel('z'); qt1 = title('Magnetization (M)');
    [a,b] = view();
    view(0,90);

subplot(122);
    q2 = quiver3(X,Y,Z, Mx(:,:,:,1), My(:,:,:,1), Mz(:,:,:,1), .5);
    axis tight; zlim ([-zdim, zdim]);
    daspect([1 1 1]);
    xlabel('x'); ylabel('y'); zlabel('z'); qt2 = title('Magnetization (M)');
    view(45,10);
    %view(3)


for i = 1:1:tdim
    set(q1, 'Udata', Mx(:,:,:,i));
    set(q1, 'Vdata', My(:,:,:,i));
    set(q1, 'Wdata', Mz(:,:,:,i));
    set(q2, 'Udata', Mx(:,:,:,i));
    set(q2, 'Vdata', My(:,:,:,i));
    set(q2, 'Wdata', Mz(:,:,:,i));
    set(qt1, 'string', ['M(t = ', num2str(time(i)), ')']);
    set(qt2, 'string', ['M(t = ', num2str(time(i)), ')']);
    drawnow;
end


%avgMx = squeeze(sum(sum(Mx)));
%avgMy = squeeze(sum(sum(My)));
%avgMz = squeeze(sum(sum(Mz)));
%avgM = sqrt(avgMx.^2 + avgMy.^2 + avgMz.^2);

%size(avgMx)

%subplot(222);
    %plot(time, avgMx, time, avgMy, time, avgMz, time, avgM);
    %legend('Mx', 'My', 'Mz', 'M');
    %xlabel('time'); title('Magnetization (A/m)');
