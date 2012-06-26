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

    xdim = 43;
    ydim = 43;
    zdim = 1;
    if(zdim == 1)       zslice = 1;
    elseif(zdim == 3)   zslice = 2;
    elseif(zdim == 4)   zslice = 3;
    elseif(zdim == 5)   zslice = 3;
    elseif(zdim == 6)   zslice = 3;
    end
    x = 0:xdim-1;
    y = 0:ydim-1;
    z = zslice:zslice;
    [X,Y,Z] = meshgrid(x,y,z);

    for bigindex = 0:0
        start_tindex = 5000*bigindex %1
        tdim = length(time);
        %start_tindex = 0;
        start_line = ydim*xdim*zdim*(start_tindex);
        required_lines = ydim*xdim*zdim*tdim;
        size(tindex)
        %tdim = 1;
    %return

    Ms = 8.6e5;
    Mfile = [folder_name, '/Mdynamics.dat'];
    %system(['tail -n +' num2str(start_line) ' Mdynamics.dat | head -n' num2str(required_lines) ' > Mdynamics1.dat']);
    %M_yxzt = load('Mdynamics1.dat');
    disp 'loading file...'; tic
    M_yxzt = load(Mfile) / Ms;
    disp 'file loaded!'; toc
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
    Mx = Mx(1:sf:end, 1:sf:end, 1:sf:end, :);
    My = My(1:sf:end, 1:sf:end, 1:sf:end, :);
    Mz = Mz(1:sf:end, 1:sf:end, 1:sf:end, :);
    animation_skip = 5;
    alfa = 1;

    fig = figure; set(fig, 'name', folder_name);
    set(gcf, 'OuterPosition', [0 0 1280 800]);
    subplot(221);
        q1 = quiver3(X,Y,Z, Mx(:,:,:,1), My(:,:,:,1), Mz(:,:,:,1), .5);
        axis tight equal; grid off;
        xlabel('x'); ylabel('y'); zlabel('z'); qt1 = title('Magnetization (M)');
        [a,b] = view();
        view(0,90);
    %subplot(222);
        %q2 = quiver3(X,Y,Z, Mx(:,:,:,1), My(:,:,:,1), Mz(:,:,:,1), .5);
        %axis tight; grid on
        %zlim ([zslice-1, zslice+1]);
        %%daspect([1 1 .5]);
        %daspect([1 1 1]);
        %xlabel('x'); ylabel('y'); zlabel('z'); qt2 = title('Magnetization (M)');
        %view(0,0);
        %set(qt1, 'string', ['M(t = ', num2str(time(start_tindex+1)), ')']);
        %set(qt2, 'string', ['M(t = ', num2str(time(start_tindex+1)), ')']);
        %%view(3)
    subplot(223);
        q3 = imagesc(X(:,1),Y(1,:), Mx(:,:,:,1).^alfa);
        xlabel('x'); ylabel('y'); title('Magnetization (Mx)');
        axis  equal tight xy;
        grid off; colorbar; colormap(jet);
        caxis([-1, 1])
    subplot(224);
        q4 = imagesc(X(:,1),Y(1,:), My(:,:,:,1).^alfa);
        xlabel('x'); ylabel('y'); title('Magnetization (My)');
        axis  equal tight xy;
        grid off; colorbar;
        caxis([-1, 1])
    subplot(222);
        q5 = imagesc(X(:,1),Y(1,:), Mz(:,:,:,1));
        xlabel('x'); ylabel('y'); title('Magnetization (Mz)');
        axis  equal tight xy;
        grid off; colorbar;
        caxis([-1, 1])


    kk = 0;
    for i = [1:animation_skip:tdim, tdim]
        set(q1, 'Udata', Mx(:,:,:,i));
        set(q1, 'Vdata', My(:,:,:,i));
        set(q1, 'Wdata', Mz(:,:,:,i));
        %set(q2, 'Udata', Mx(:,:,:,i));
        %set(q2, 'Vdata', My(:,:,:,i));
        %set(q2, 'Wdata', Mz(:,:,:,i));
        set(qt1, 'string', ['M(t = ', num2str(time(i+start_tindex)), ')']);
        %set(qt2, 'string', ['M(t = ', num2str(time(i+start_tindex)), ')']);
        set(q3, 'Cdata', Mx(:,:,:,i).^alfa);
        set(q4, 'Cdata', My(:,:,:,i).^alfa);
        set(q5, 'Cdata', Mz(:,:,:,i));
        drawnow;
        kk = kk + 1;
        fnum = sprintf('Mcolors%05d', kk)
%         print(gcf, fnum, '-depsc');
    end

%     print(gcf, ['M' num2str(start_tindex)], '-dpdf');
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

end % function
