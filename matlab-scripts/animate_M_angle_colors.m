function animate_M_angle_colors(folder_name)

clc

%% decide to post-process or not
datfilename = [folder_name, '/Mdynamics.dat'];
matfilename = [folder_name, '/Mdynamics-pp.mat'];
dateDAT = dir(datfilename);
dateMAT = dir(matfilename);
if length(dateMAT) && (dateMAT.datenum > dateDAT.datenum)
    disp([matfilename,' exists already'])
    doReadFile = false;
    disp(['loading ',matfilename,' ...'])
    load(matfilename)
else
    disp([matfilename,' is old or does not exist'])
    doReadFile = true;
end

%% Post-processing dat file
if doReadFile
    disp(['post-processing ', datfilename,' ...'])
    %% load dynamics
    filename = [folder_name, '/dynamics.dat'];
    dynamics = load(filename);
    tindex = dynamics(:,1);
    time = dynamics(:,2);
    % dt = dynamics(:,3);
    % E = dynamics(:,4);
    % Mx = dynamics(:,5);
    % My = dynamics(:,6);
    % Mz = dynamics(:,7);
    % M  = dynamics(:,8);
    % torque  = dynamics(:,9);
    clear dynamics
    tdim = length(time)

    %% read Mdynamics
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
    Mx = zeros(Ny,Nx,tdim);
    My = zeros(Ny,Nx,tdim);
    Mz = zeros(Ny,Nx,tdim);
    disp(['rearranging data ...']);
    for i = 0:tdim-1
        m = C(1, i*Nx*Ny+1:(i+1)*Nx*Ny);
        Mx(:,:,i+1) = reshape(m,Nx,Ny)';
        m = C(2, i*Nx*Ny+1:(i+1)*Nx*Ny);
        My(:,:,i+1) = reshape(m,Nx,Ny)';
        m = C(3, i*Nx*Ny+1:(i+1)*Nx*Ny);
        Mz(:,:,i+1) = reshape(m,Nx,Ny)';
    end
    %% reduce variables
    Ms = 8.6e5;
    mask = ((Mx(:,:,1).^2+My(:,:,1).^2+Mz(:,:,1).^2) ~= 0);
    theta = acos(Mz/Ms) * 180/pi;
    phi = atan2(My,Mx) * 180/pi;
    clear Mx My Mz
    clear C num fid i m dateDAT dateMAT filename doReadFile
    whos
    disp(['saving ', matfilename,' ...']);
    save(matfilename)
end


%% Visualization params
doAnimation = true;
animation_skip = 2;
alfa = 1;
sf = 1;


%% subsample
x = x(1:sf:end);
y = y(1:sf:end);
[X,Y] = meshgrid(x,y);
Z = 0*X;
Nx = length(x);
Ny = length(y);
theta = theta(1:sf:end, 1:sf:end, :);
phi   = phi  (1:sf:end, 1:sf:end, :);
mask  = mask (1:sf:end, 1:sf:end, :);

%% animation setup
Mx = repmat(mask,[1 1 tdim]) .* sin(theta*pi/180) .* cos(phi*pi/180);
My = repmat(mask,[1 1 tdim]) .* sin(theta*pi/180) .* sin(phi*pi/180);
Mz = repmat(mask,[1 1 tdim]) .* cos(theta*pi/180);

clf % fig = figure; 
% set(gcf, 'OuterPosition', [0 0 1280 800]);
subplot(211);
    q1 = quiver3(X,Y,Z+.25, Mx(:,:,1), My(:,:,1), Mz(:,:,1), .5);
    set(q1,'color','black','linewidth',0.5);
    grid off;
    xlabel('x'); ylabel('y'); title('Magnetization (\phi - angle in plane)');
    view(3);
%     set(gca,'visible','off')
    hold on
    q3 = pcolor(X(1,:),Y(:,1), phi(:,:,1).^alfa);
    hold off;
    shading interp;
    axis equal
    axis  equal tight xy;
    %     colorbar; 
    grid off; 
    colormap(hsv);
    caxis([-180,180])
    
h2 = subplot(212);
    q2 = plot(time(1:2)/1e-9,squeeze(Mx(round(Ny/2),round(Nx/2),1:2)));
    hold on
    q2b = plot(time(1:2)/1e-9,squeeze(Mx(round(Ny/4),round(Nx/2),1:2)),'r');
    hold off;
    set(gca,'ylim',[-1 1])
    legend('Mx @ center', 'Mx @ 1/4 of diameter')
    xlabel('time [ns]');
    qt2 = title('Magnetiztion');

%     return
    
%% animation loop
    if ~doAnimation
        animation_skip = tdim-1;
    end    
    kk = 0;
    for i = [2:animation_skip:tdim, tdim]
        set(q1, 'Udata', Mx(:,:,i));
        set(q1, 'Vdata', My(:,:,i));
        set(q1, 'Wdata', Mz(:,:,i));

        set(q3, 'Xdata', X(1,:));
        set(q3, 'Ydata', Y(:,1));
        set(q3, 'Cdata', phi(:,:,i).^alfa);

        set(q2, 'xdata', time(1:i)/1e-9);
        set(q2, 'ydata', squeeze(Mx(round(Ny/2),round(Nx/2),1:i)));
        set(q2b, 'xdata', time(1:i)/1e-9);
        set(q2b, 'ydata', squeeze(Mx(round(Ny/4),round(Nx/2),1:i)));
        set(qt2, 'string', ['M(t = ', num2str(time(i+0)), ')']);
        set(h2,'xlim',[time(1) time(i)]/1e-9);
        
        if doAnimation        
            drawnow;
        end
        
        kk = kk + 1;
        fnum = sprintf('Mcolors%05d', kk);
%         disp(fnum)
%         print(gcf, fnum, '-depsc');

%         if kk == 1
%             break
%         end
    end
    print(gcf, 'Manimate', '-dpdf');

% keyboard

    