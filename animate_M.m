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
clear dynamics

for bigindex = 0:5
    animation_skip = 1;
    start_tindex = 5000*bigindex %1
    tdim = 1 %length(time);
    zdim = 5;
    ydim = 32;
    xdim = ydim;
    zslice = 3;
    x = 0:xdim-1;
    y = 0:ydim-1;
    z = zslice:zslice;
    [X,Y,Z] = meshgrid(x,y,z);
    start_line = ydim*xdim*zdim*(start_tindex);
    required_lines = ydim*xdim*zdim*tdim;
%return

%system(['tail -' num2str(required_lines) ' Mdynamics.dat > Mdynamics1.dat']);
system(['tail -n +' num2str(start_line) ' Mdynamics.dat | head -n' num2str(required_lines) ' > Mdynamics1.dat']);
M_yxzt = load('Mdynamics1.dat');
    M_yxzt = M_yxzt(1:required_lines,:);
    M_yxzt = reshape(M_yxzt', 3,ydim,xdim,zdim, tdim);
    Mx = shiftdim(M_yxzt(1,:,:,zslice,:), 1);
    My = shiftdim(M_yxzt(2,:,:,zslice,:), 1);
    Mz = shiftdim(M_yxzt(3,:,:,zslice,:), 1);
clear M_yxzt

%figure;
%set(gcf, 'OuterPosition', [0 0 1280 800]);
subplot(121);
    q1 = quiver3(X,Y,Z, Mx(:,:,:,1), My(:,:,:,1), Mz(:,:,:,1), .5);
    axis tight equal; grid off;
    xlabel('x'); ylabel('y'); zlabel('z'); qt1 = title('Magnetization (M)');
    [a,b] = view();
    view(0,90);

subplot(122);
    q2 = quiver3(X,Y,Z, Mx(:,:,:,1), My(:,:,:,1), Mz(:,:,:,1), .5);
    axis tight; grid on
    zlim ([zslice-1, zslice+1]);
    daspect([1 1 .5]);
    xlabel('x'); ylabel('y'); zlabel('z'); qt2 = title('Magnetization (M)');
    view(45,10);
    set(qt1, 'string', ['M(t = ', num2str(time(start_tindex+1)), ')']);
    set(qt2, 'string', ['M(t = ', num2str(time(start_tindex+1)), ')']);
    %view(3)

for i = 1:animation_skip:tdim
    set(q1, 'Udata', Mx(:,:,:,i));
    set(q1, 'Vdata', My(:,:,:,i));
    set(q1, 'Wdata', Mz(:,:,:,i));
    set(q2, 'Udata', Mx(:,:,:,i));
    set(q2, 'Vdata', My(:,:,:,i));
    set(q2, 'Wdata', Mz(:,:,:,i));
    set(qt1, 'string', ['M(t = ', num2str(time(i+start_tindex)), ')']);
    set(qt2, 'string', ['M(t = ', num2str(time(i+start_tindex)), ')']);
    drawnow;
end

print(gcf, ['M' num2str(start_tindex)], '-dpdf');
end % bigindex



%avgMx = squeeze(sum(sum(Mx)));
%avgMy = squeeze(sum(sum(My)));
%avgMz = squeeze(sum(sum(Mz)));
%avgM = sqrt(avgMx.^2 + avgMy.^2 + avgMz.^2);

%size(avgMx)

%subplot(222);
    %plot(time, avgMx, time, avgMy, time, avgMz, time, avgM);
    %legend('Mx', 'My', 'Mz', 'M');
    %xlabel('time'); title('Magnetization (A/m)');
