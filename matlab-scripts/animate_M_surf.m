function animate_M_surf(folder_name)

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
%         Hext  = dynamics(:,10);
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

    for bigindex = 0:0
        start_tindex = 5000*bigindex %1
        tdim = length(time);
        %start_tindex = 0;
        start_line = ydim*xdim*zdim*(start_tindex);
        required_lines = ydim*xdim*zdim*tdim;
        size(tindex)
        %tdim = 1;
    %return

    Mfile = [folder_name, '/Mdynamics.dat'];
    %system(['tail -n +' num2str(start_line) ' Mdynamics.dat | head -n' num2str(required_lines) ' > Mdynamics1.dat']);
    %M_yxzt = load('Mdynamics1.dat');
    disp 'loading file...'; tic
    M_yxzt = load(Mfile);
    disp 'file loaded!'; toc
        M_yxzt = M_yxzt(1:required_lines,:);
        M_yxzt = reshape(M_yxzt', 3,ydim,xdim,zdim, tdim);
        Mx = shiftdim(M_yxzt(1,:,:,zslice,:), 1);
        My = shiftdim(M_yxzt(2,:,:,zslice,:), 1);
        Mz = shiftdim(M_yxzt(3,:,:,zslice,:), 1);
    clear M_yxzt

% subsample
    sf = 2;
    X = X(1:sf:end, 1:sf:end, 1:sf:end);
    Y = Y(1:sf:end, 1:sf:end, 1:sf:end);
    Z = Z(1:sf:end, 1:sf:end, 1:sf:end);
    Mx = Mx(1:sf:end, 1:sf:end, 1:sf:end, :);
    My = My(1:sf:end, 1:sf:end, 1:sf:end, :);
    Mz = Mz(1:sf:end, 1:sf:end, 1:sf:end, :);
    animation_skip = 20;

    fig = figure; set(fig, 'name', folder_name);
    set(gcf, 'OuterPosition', [0 0 1280 800]);

    %subplot(221);
        %q1 = quiver3(X,Y,Z, Mx(:,:,:,1), My(:,:,:,1), Mz(:,:,:,1), .5);
        %axis tight equal; grid off;
        %xlabel('x'); ylabel('y'); zlabel('z'); qt1 = title('Magnetization (M)');
        %[a,b] = view();
        %view(0,90);
    subplot(221);
%         q1 = plot(time(1), Hext(1));
%         grid on;
%         xlabel('time'); ylabel('H_{ext}'); qt1 = title('External Field (Hext)');
%         axis([time(1), time(end) -5e5 5e5])
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
        q3 = surf(X, Y, Mx(:,:,:,1)); shading faceted
        %[C, q3] = contourf(X, Y, Mx(:,:,:,1)); %shading interp
        xlabel('x'); ylabel('y'); title('Magnetization (Mx)');
        grid off; colorbar; colormap(jet);
        %axis square
        zlim(8.6e5*[-1, 1])
        caxis(8.6e5*[-1, 1])
        axis square
        %view(0,0)
    subplot(224);
        q4 = surf(X, Y, My(:,:,:,1)); shading faceted
        %[C, q4] = contourf(X, Y, My(:,:,:,1)); %shading interp
        xlabel('x'); ylabel('y'); title('Magnetization (My)');
        grid off; colorbar;
        %axis square
        zlim(8.6e5*[-1, 1])
        caxis(8.6e5*[-1, 1])
        axis square
        %view(0,0)
    subplot(222);
        q5 = surf(X, Y, Mz(:,:,:,1)); shading faceted
        %[C, q5] = contourf(X, Y, Mz(:,:,:,1)); %shading interp
        xlabel('x'); ylabel('y'); title('Magnetization (Mz)');
        grid off; colorbar;
        %axis square
        zlim(8.6e5*[-1, 1])
        caxis(8.6e5*[-1, 1])
        axis square
    %subplot(224);
        %q4 = imagesc(X(:,1),Y(1,:), My(:,:,:,1));
        %xlabel('x'); ylabel('y'); title('Magnetization (My)');
        %axis  equal tight xy;
        %grid off; colorbar;
        %caxis(8.6e5*[-1, 1])

    %keyboard

    %keyboard
    pause(2)
    kk = 0;
    for i = [1:animation_skip:tdim, tdim]
        %set(q1, 'Udata', Mx(:,:,:,i));
        %set(q1, 'Vdata', My(:,:,:,i));
        %set(q1, 'Wdata', Mz(:,:,:,i));
%         set(q1, 'Xdata', time(1:i));
%         set(q1, 'Ydata', Hext(1:i));
        %set(q2, 'Udata', Mx(:,:,:,i));
        %set(q2, 'Vdata', My(:,:,:,i));
        %set(q2, 'Wdata', Mz(:,:,:,i));
%         set(qt1, 'string', ['M(t = ', num2str(time(i+start_tindex)), ')']);
        %set(qt2, 'string', ['M(t = ', num2str(time(i+start_tindex)), ')']);
        set(q3, 'Zdata', Mx(:,:,:,i));
        set(q4, 'Zdata', My(:,:,:,i));
        set(q5, 'Zdata', Mz(:,:,:,i));
        %set(q5, 'Cdata', Mz(:,:,:,i));
        drawnow;
        kk = kk + 1;
        fnum = sprintf('Mcont%05d', kk)
        print(gcf, fnum, '-depsc');
        %print(gcf, ['M' num], '-depsc');
    end

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


    %keyboard
end % function
