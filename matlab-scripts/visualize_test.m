function animate_M(folder_name)

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

    %fig = figure; set(fig, 'name', folder_name);
    set(gcf, 'OuterPosition', [1 1 1280 800]);
    subplot(221);
        plot(time, E, '-');
        grid on;
        xlabel('time'); title('Energy (eV)');

    subplot(222);
        plot(time, Mx, time, My, time, Mz, time, M);
%         legend('Mx', 'My', 'Mz', 'M');
%         plot(time, Mx, time, My);
%         legend('Mx', 'My');
        ylim([-1 1]);
%         ylim([-0.02 0.02]);
        grid on;
        xlabel('time'); title('Magnetization (A/m)');



    subplot(223);
        %plot(time, torque, '-');
        semilogy(time, torque, '-');
        grid on;
        xlabel('time'); title('Normalized maximum Torque M \times H / Ms^2');

    subplot(223);
%         plot(time, Hext, '.-');
%         grid on;
%         xlabel('time'); title('Hext (A/m)');

   subplot(224);
        %plot(time, dt, '-');
        semilogy(time, dt, '-');
        %axis tight;
        grid on;
        xlabel('time'); title('Time step');

    print(gcf, ['dynamics'], '-depsc');

end % function
