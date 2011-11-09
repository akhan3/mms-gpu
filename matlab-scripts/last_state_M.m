function last_state_M(folder_name)

datfilename = [folder_name, '/Mlatest.dat'];
doReadFile = true;

%% Post-processing dat file
if doReadFile
    disp(['post-processing ', datfilename,' ...'])
    %% load dynamics
    filename = [folder_name, '/dynamics.dat'];
    dynamics = load(filename);
    tindex = dynamics(:,1);
    time = dynamics(:,2);
    dt = dynamics(:,3);
    E = dynamics(:,4);
    Mx_avg = dynamics(:,5);
    My_avg = dynamics(:,6);
    Mz_avg = dynamics(:,7);
    M_avg  = dynamics(:,8);
    torque  = dynamics(:,9);
    clear dynamics
    tdim = length(time)
    %% read Mlatest
    fid = fopen(datfilename);
    [C,num] = textscan(fid, 'Nx=%d\nNy=%d\n');
    Nx = cell2mat(C(1));
    Ny = cell2mat(C(2));
    x = 0:Nx-1;
    y = 0:Ny-1;
    [X,Y] = meshgrid(x,y);
    Z = 0*X;
    [C,num] = textscan(fid, '# start field');
    [C,num] = textscan(fid, '%f');
    disp(['read from ', datfilename,' done'])
    C = cell2mat(C);
    C = reshape(C,3,length(C)/3);
    disp(['rearranging data ...']);
    Mx(:,:) = reshape(C(1,:),Nx,Ny)';
    My(:,:) = reshape(C(2,:),Nx,Ny)';
    Mz(:,:) = reshape(C(3,:),Nx,Ny)';
    Ms = 8.6e5;
end


%% subsample
sf = 2;
x = x(1:sf:end);
y = y(1:sf:end);
[X,Y] = meshgrid(x,y);
Z = 0*X;
Nx = length(x);
Ny = length(y);
Mx = Mx(1:sf:end, 1:sf:end);
My = My(1:sf:end, 1:sf:end);
Mz = Mz(1:sf:end, 1:sf:end);

%% plot
clf % fig = figure; set(fig, 'name', filename);
set(gcf, 'OuterPosition', [0 0 1280 800]);

subplot(221);
    q1 = quiver3(X,Y,Z, Mx(:,:,1), My(:,:,1), Mz(:,:,1), .5);
%     q1 = quiver3(X,Y,Z, Mx(:,:,1), My(:,:,1), zeros(size(Mz(:,:,1))), .5);
    set(q1,'color','black','linewidth',0.2);
    grid off;
    xlabel('x'); ylabel('y'); title('Magnetization (\phi - angle in plane)');
    view(3);
    set(gca,'visible','off')
    hold on
    phi = atan2(My,Mx) * 180/pi;
    q3 = pcolor(X(1,:),Y(:,1), phi(:,:,1));
    hold off;
    shading interp;
    axis equal
    axis  equal tight xy;
    grid off;
    colormap(hsv);
    % colorbar
    caxis([-180,180])


subplot(222);
    plot(time/1e-9, E, '-');
    grid on;
    xlabel('time [ns]'); title('Energy (eV)');
    

subplot(223);
    plot(time/1e-9, Mx_avg, time/1e-9, My_avg, time/1e-9, Mz_avg, time/1e-9, M_avg);
    %     legend('Mx', 'My', 'Mz', 'M');
    ylim([-1 1]);
    grid on;
    xlabel('time [ns]'); title('Magnetization (A/m)');


subplot(224);
    plot(time/1e-9, torque, '-');
%     semilogy(time/1e-9, torque, '-');
    grid on;
    xlabel('time [ns]'); title('Normalized maximum Torque M \times H / Ms^2');


print(gcf, ['Mlatest'], '-dpdf');
return


