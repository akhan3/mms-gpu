 %initial loading and setup
    clear

load dynamics.dat
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

%set(gcf, 'OuterPosition', [0 0 1280 800]);
subplot(221);
    plot(time, E, '.-');
    %axis tight;
    xlabel('time'); title('Energy (eV)');

subplot(222);
    plot(time, Mx, time, My, time, Mz, time, M);
    %axis tight;
    ylim([-1 1]);
    legend('Mx', 'My', 'Mz', 'M');
    xlabel('time'); title('Magnetization (A/m)');



subplot(223);
    %plot(time, torque, '-');
    semilogy(time, torque, '.-');
    %axis tight;
    xlabel('time'); title('Normalized maximum Torque M \times H / Ms^2');

subplot(224);
    %plot(time, dt, '-');
    semilogy(time, dt, '.-');
    %axis tight;
    xlabel('time'); title('Time step');
