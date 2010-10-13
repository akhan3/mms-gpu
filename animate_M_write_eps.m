function animate_M(folder_name)

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

    xdim = 42;
    ydim = 26;
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

    tdim = length(time);

    Mfile = [folder_name, '/Mdynamics.dat'];
    M_yxzt = load(Mfile);
        M_yxzt = M_yxzt(1:end, :);
        M_yxzt = reshape(M_yxzt', 3,ydim,xdim,zdim, tdim);
        Mx = shiftdim(M_yxzt(1,:,:,zslice, :), 1);
        My = shiftdim(M_yxzt(2,:,:,zslice, :), 1);
        Mz = shiftdim(M_yxzt(3,:,:,zslice, :), 1);
    clear M_yxzt

% subsample
    sf = 1;
    X = X(1:sf:end, 1:sf:end, 1:sf:end);
    Y = Y(1:sf:end, 1:sf:end, 1:sf:end);
    Z = Z(1:sf:end, 1:sf:end, 1:sf:end);
    Mx = Mx(1:sf:end, 1:sf:end, 1:sf:end, :);
    My = My(1:sf:end, 1:sf:end, 1:sf:end, :);
    Mz = Mz(1:sf:end, 1:sf:end, 1:sf:end, :);
    animation_skip = 10;

    fig = figure; set(fig, 'name', folder_name);
    set(gcf, 'OuterPosition', [0 0 1280 800]);
    q1 = quiver3(X,Y,Z, Mx(:,:,:,1), My(:,:,:,1), Mz(:,:,:,1), .5);
    set(gca, 'visible', 'off');
    axis tight equal; grid off;
    view(0,90);
    frame_count = 1;

    %tindex_start = 1;
    %tindex_end = 335;
    %start_line = ydim*xdim*zdim*(start_tindex);
    %required_lines = ydim*xdim*zdim*tdim;

    %for i = [ 1:1:335,  335:10:tdim,  tdim ]
        %set(q1, 'Udata', Mx(:,:,:,i));
        %set(q1, 'Vdata', My(:,:,:,i));
        %set(q1, 'Wdata', Mz(:,:,:,i));
        %drawnow;
        %print(gcf, [folder_name, '/M' num2str(frame_count)], '-depsc');
        %%print(gcf, [folder_name, '/M' num2str(frame_count)], '-dpng');
        %fprintf('frame_count=%d, i=%d\n', frame_count, i);
        %frame_count = frame_count+1;
    %end

    torque_max = 1e1;
    skip_min = 1;
    torque_min = 1e-3;
    skip_max = 100;
    skip_slope = (skip_max - skip_min) / (log10(torque_min) - log10(torque_max));

    i = 1;
    while i <= tdim
        if(torque(i) > torque_max)
            skip = skip_min;
        else
            skip = floor(skip_slope * (log10(torque(i)) - log10(torque_max)) + skip_min);
        end

        %elseif (torque(i) > 0.01)   skip = 5;
        %elseif (torque(i) > 0.001)  skip = 10;
        %else                        skip = 10;
        %end

        set(q1, 'Udata', Mx(:,:,:,i));
        set(q1, 'Vdata', My(:,:,:,i));
        set(q1, 'Wdata', Mz(:,:,:,i));
        drawnow;
        str_frame_count = sprintf('%04d', frame_count);
        print(gcf, [folder_name, '/M' str_frame_count], '-depsc');
        fprintf('frame_count=%d, i=%d, skip=%d, torque=%f\n', frame_count, i, skip, torque(i));
        frame_count = frame_count + 1;
        i = i + skip;
    end

end
