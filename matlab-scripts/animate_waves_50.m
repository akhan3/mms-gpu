function animate_M(folder_name)
if (1)
%     folder_name = '../results/sim_spin_50x50x3'

    dynamics = load([folder_name, '/dynamics.dat']);
    dynamics = dynamics(1:end-1,:);
        tindex = dynamics(:,1);
        time = dynamics(:,2);
        dt = dynamics(:,3);
        E = dynamics(:,4);
        Mx = dynamics(:,5);
        My = dynamics(:,6);
        Mz = dynamics(:,7);
        M  = dynamics(:,8);
        torque  = dynamics(:,9);
    clear dynamicss

    Ms = 8.6e5;
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

    start_tindex = 0;
    tdim = length(time);
    required_lines = ydim*xdim*zdim*tdim;
    size(tindex)

    tic
        Mfile = [folder_name, '/Mdynamics.dat'];
        M_yxzt = load(Mfile);
            M_yxzt = M_yxzt(1:required_lines,:) ./ Ms;
            M_yxzt = reshape(M_yxzt', 3,ydim,xdim,zdim, tdim);
            Mx = shiftdim(M_yxzt(1,:,:,zslice,:), 1);
            My = shiftdim(M_yxzt(2,:,:,zslice,:), 1);
            Mz = shiftdim(M_yxzt(3,:,:,zslice,:), 1);
        clear M_yxzt
    toc

    dia = xdim * 5e-9;
    xbar = (linspace(-dia/2, dia/2, xdim));

end





% subsample
    animation_skip = 10;
    sf = 1;
    X = X(1:sf:end, 1:sf:end, 1:sf:end);
    Y = Y(1:sf:end, 1:sf:end, 1:sf:end);
    Z = Z(1:sf:end, 1:sf:end, 1:sf:end);
    Mx = Mx(1:sf:end, 1:sf:end, 1:sf:end, :);
    My = My(1:sf:end, 1:sf:end, 1:sf:end, :);
    Mz = Mz(1:sf:end, 1:sf:end, 1:sf:end, :);
    power = 1;

    fig = figure; set(fig, 'name', folder_name);
    set(gcf, 'OuterPosition', [1 1 1280 800]);

%     subplot(221);
%         q1 = plot(time(1), Hext(1));
%         grid on;
%         xlabel('time (s)'); ylabel('H_{ext}'); qt1 = title('External Field Hext(t = 0)');
%         axis([time(1), time(end) -5e5 5e5])

    subplot(222);
        Mx_avgline = Mx(end/2,:,:,1).^2;
        q2 = semilogy(xbar/1e-9, Mx_avgline);
        xlabel('position (nm)'); ylabel('M_x^2 / M_s^2'); title('Running average of squared Magnetization (Mx)');
        ylim([1e-6 1]);
        grid on;

    subplot(223);
        q3 = plot(xbar/1e-9, Mx(end/2,:,:,1).^power)
        xlabel('position (nm)'); ylabel('M_x / M_s'); title('Magnetization (Mx)');
        ylim([-1 1]);
        grid on;

    subplot(224);
        q4 = plot(xbar/1e-9, My(end/2,:,:,1).^power)
        xlabel('position (nm)'); ylabel('M_y / M_s'); title('Magnetization (My)');
        ylim([-1 1]);
        grid on;

    pause(1)
    %return

    kk = 0;
    %for i = [1:animation_skip:tdim, tdim]
    for i = [1:tdim, tdim]
        time_str = sprintf('%1.2e', time(i+start_tindex));
        Mx_avgline = ((i-1) * Mx_avgline + Mx(end/2,:,:,i).^2) / i;
%         set(q1, 'Xdata', time(1:i));
%         set(q1, 'Ydata', Hext(1:i));
%         set(qt1, 'string', ['External Field Hext(t = ', time_str, ')']);
        %set(q2, 'Ydata', Mz(end/2,:,:,i).^power);
        set(q2, 'Ydata', Mx_avgline/max(Mx_avgline));
        set(q3, 'Ydata', Mx(end/2,:,:,i).^power);
        set(q4, 'Ydata', My(end/2,:,:,i).^power);
        if(mod(i, animation_skip) == 0)
            kk = kk + 1;
            num = sprintf('%05d', kk);
            drawnow; disp(num)
            print(gcf, ['M_waves_' num], '-depsc');
        end
    end

    set(q2, 'Ydata', Mx_avgline/max(Mx_avgline));
    drawnow; disp(num)




    %avgMx = squeeze(sum(sum(Mx)));
    %avgMy = squeeze(sum(sum(My)));
    %avgMz = squeeze(sum(sum(Mz)));
    %avgM = sqrt(avgMx.^2 + avgMy.^2 + avgMz.^2);

    %size(avgMx)

    %subplot(222);
        %plot(time, avgMx, time, avgMy, time, avgMz, time, avgM);
        %legend('Mx', 'My', 'Mz', 'M');
        %xlabel('time'); title('Magnetization (A/m)');


%end % function
