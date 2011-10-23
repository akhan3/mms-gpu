function visualize_test(folder_name)

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
        %Hext  = dynamics(:,10);
    clear dynamics
    Npoints = length(time);
    fprintf('%d points in %g ns\n', Npoints, time(end)/1e-9);

%     fig = figure; set(fig, 'name', folder_name);
    set(gcf, 'OuterPosition', [1 1 1280 800]);
    subplot(221);
        plot(time/1e-9, E, '-');
        grid on;
        xlabel('time [ns]'); title('Energy (eV)');

    subplot(222);
        plot(time/1e-9, Mx, time/1e-9, My, time/1e-9, Mz, time/1e-9, M);
%         legend('Mx', 'My', 'Mz', 'M');
        ylim([-1 1]);
        grid on;
        xlabel('time [ns]'); title('Magnetization (A/m)');



    subplot(223);
        %plot(time/1e-9, torque, '-');
        semilogy(time/1e-9, torque, '-');
        grid on;
        xlabel('time [ns]'); title('Normalized maximum Torque M \times H / Ms^2');

   subplot(224);
        %plot(time/1e-9, dt, '-');
        semilogy(time/1e-9, dt, '-');
        %axis tight;
        grid on;
        xlabel('time [ns]'); title('Time step');

    print(gcf, ['dynamics'], '-depsc');

end % function
