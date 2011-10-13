function last_state_M(folder_name)

    dynamics = load([folder_name, '/dynamics.dat']);
        tindex = dynamics(:,1);
        time = dynamics(:,2);
        dt = dynamics(:,3);
        E = dynamics(:,4);
        Mx = dynamics(:,5);
        My = dynamics(:,6);
        Mz = dynamics(:,7);
        M  = dynamics(:,8);
        torque  = dynamics(:,9);
    clear dynamics

    xdim = 50;
    ydim = 50;
    zdim = 3;
    if(zdim == 3)       zslice = 2;
    elseif(zdim == 4)   zslice = 3;
    elseif(zdim == 5)   zslice = 3;
    elseif(zdim == 6)   zslice = 3;
    end
    x = 0:xdim-1;
    y = 0:ydim-1;
    z = zslice:zslice;
    [X,Y,Z] = meshgrid(x,y,z);

    ti = tindex(end)-1;
    tdim = 1;
    start_line = ti * ydim*xdim*zdim;
    required_lines = ydim*xdim*zdim*tdim;

    Mfile = [folder_name, '/Mdynamics.dat'];
    system(['tail -n +' num2str(start_line) ' ' Mfile ' | head -n' num2str(required_lines) ' > Mlast.dat']);
    M_yxzt = load('Mlast.dat');
        M_yxzt = M_yxzt(1:required_lines,:);
        M_yxzt = reshape(M_yxzt', 3,ydim,xdim,zdim, tdim);
        Mx = shiftdim(M_yxzt(1,:,:,zslice,:), 1);
        My = shiftdim(M_yxzt(2,:,:,zslice,:), 1);
        Mz = shiftdim(M_yxzt(3,:,:,zslice,:), 1);
    clear M_yxzt

% subsample
    sf = 1;
    X = X(1:sf:end, 1:sf:end, 1:sf:end);
    Y = Y(1:sf:end, 1:sf:end, 1:sf:end);
    Z = Z(1:sf:end, 1:sf:end, 1:sf:end);
    Mx = Mx(1:sf:end, 1:sf:end, 1:sf:end, 1);
    My = My(1:sf:end, 1:sf:end, 1:sf:end, 1);
    Mz = Mz(1:sf:end, 1:sf:end, 1:sf:end, 1);

    fig = figure; set(fig, 'name', folder_name);
    %fig = figure; set(fig, 'name', [num2str(xdim), 'x', num2str(ydim), 'x', num2str(zdim)]);
    set(gcf, 'OuterPosition', [0 0 1280 800]);
    q1 = quiver3(X,Y,Z, Mx(:,:,:,1), My(:,:,:,1), Mz(:,:,:,1), .5);
    set(gca, 'visible', 'off');
    %set(gca, 'visible', 'off');
    axis tight equal; grid off;
    xlabel('x'); ylabel('y'); zlabel('z'); qt1 = title('Magnetization (M)');
    [a,b] = view();
    view(0,90);


    print(gcf, ['M', 'state'], '-depsc');
    return








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
        %daspect([1 1 .5]);
        daspect([1 1 1]);
        xlabel('x'); ylabel('y'); zlabel('z'); qt2 = title('Magnetization (M)');
        view(0,0);
        set(qt1, 'string', ['M(t = ', num2str(time(end)), ') E = ', num2str(E(end)), 'eV']);
        set(qt2, 'string', ['M(t = ', num2str(time(end)), ') E = ', num2str(E(end)), 'eV']);
        %view(3)


end % function
