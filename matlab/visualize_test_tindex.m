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

    fig = figure; set(fig, 'name', folder_name);
    set(gcf, 'OuterPosition', [0 0 1280 800]);
    subplot(221);
        plot(tindex, E, '.-');
        %axis tight;
        grid on;
        xlabel('tindex'); title('Energy (eV)');

    subplot(222);
        plot(tindex, Mx, tindex, My, tindex, Mz, tindex, M);
        %axis tight;
        ylim([-1 1]);
        grid on;
        %legend('Mx', 'My', 'Mz', 'M');
        xlabel('tindex'); title('Magnetization (A/m)');



    subplot(223);
        %plot(tindex, torque, '-');
        semilogy(tindex, torque, '.-');
        %axis tight;
        grid on;
        xlabel('tindex'); title('Normalized maximum Torque M \times H / Ms^2');

    subplot(224);
        %plot(tindex, dt, '-');
        semilogy(tindex, dt, '.-');
        %axis tight;
        grid on;
        xlabel('tindex'); title('Time step');

    print(gcf, ['dynamics'], '-depsc');

    %keyboard
end % function
