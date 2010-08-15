 %initial loading and setup
    clear

load dynamics.dat
    tindex = dynamics(:,1);
    time = dynamics(:,2);
    dt = dynamics(:,3);
    energy = dynamics(:,4);
    dE = dynamics(:,5); dE(1) = dE(2);
    Mx = dynamics(:,6);
    My = dynamics(:,7);
    Mz = dynamics(:,8);
    M  = dynamics(:,9);
    torque_max  = dynamics(:,10);
clear dynamics

%set(gcf, 'OuterPosition', [0 0 1280 800]);
subplot(221);
    plot(time, energy, '-');
    %axis tight;
    xlabel('time'); title('Energy (eV)');

subplot(222);
    plot(time, Mx, time, My, time, Mz, time, M);
    %axis tight;
    legend('Mx', 'My', 'Mz', 'M');
    xlabel('time'); title('Magnetization (A/m)');



subplot(223);
    %plot(time, torque_max, '-');
    semilogy(time, torque_max, '-');
    %axis tight;
    xlabel('time'); title('Normalized maximum Torque M \times H / Ms^2');

subplot(224);
    %plot(time, dt, '-');
    semilogy(time, dt, '-');
    %axis tight;
    xlabel('time'); title('Time step');


return
% load and post-process M data
    zdim = 32;
    ydim = 64;
    xdim = ydim;
    x = 0:xdim-1;
    y = 0:ydim-1;
    z = 0:zdim-1;
    [X,Y,Z] = meshgrid(x,y,z);
    Mall = load('M.dat');
    Mall = reshape(Mall', 3,ydim,xdim,zdim);
    Mx = shiftdim(Mall(1,:,:,:), 1);
    My = shiftdim(Mall(2,:,:,:), 1);
    Mz = shiftdim(Mall(3,:,:,:), 1);
    M = sqrt(Mx.^2 + My.^2 + Mz.^2);
    Ms = 8.6e5;
subplot(224);
    quiver3(X,Y,Z, Mx,My,Mz, .5);
    axis tight equal;
    xlabel('x'); ylabel('y'); zlabel('z'); title('Magnetization (M)');
    [a,b] = view();
    view(0,90);

%subplot(224);
    %quiver3(X,Y,Z, Mx,My,Mz, .5);
    %axis tight; daspect([1 1 .2]); zlim ([-zdim, zdim]);
    %xlabel('x'); ylabel('y'); zlabel('z'); title('Magnetization (M)');
    %view(0,0);

    %subplot(234); imagesc(Mx); axis image xy; caxis([-Ms Ms]); xlabel('x'); ylabel('y'); zlabel('z'); title('Mx');
    %subplot(235); imagesc(My); axis image xy; caxis([-Ms Ms]); xlabel('x'); ylabel('y'); zlabel('z'); title('My');
    %subplot(236); imagesc(Mz); axis image xy; caxis([-Ms Ms]); xlabel('x'); ylabel('y'); zlabel('z'); title('Mz');

return





%return

%load fmm_data_64x64.mat
 %draw H field in slices
    %figure('Colormap',map)
%subplot(223)
%clf
    slice_num = ceil(zdim/2);
    contourslice(H,[],[],[1, slice_num, zdim])
    axis image xy
    view(3)
    xlabel('x'); ylabel('y'); zlabel('z');
    %daspect([1 1 .02])

return

    imagesc(M(:,:,zdim/2));
    axis image xy; colorbar; hold on;
    streamslice(Hx(:,:,zdim/2), Hy(:,:,zdim/2));
    %quiver(x,y, Hx(:,:,zdim/2), Hy(:,:,zdim/2), .5)
    hold off;



%subplot(223);
    %quiver3(X,Y,Z, Mx,My,Mz, .5);
    %axis tight equal;
    %xlabel('x'); ylabel('y'); zlabel('z');
    %title('Magnetization (M)');
    %[a,b] = view();
    %view(a+90-30,b);
    %%view(0,90);

%subplot(224);
    %quiver3(X,Y,Z, Mx,My,Mz, .5);
    %axis tight equal;
    %xlabel('x'); ylabel('y'); zlabel('z');
    %title('Magnetization (M)');
    %view(90,0);

    %zz = 2;
    %imagesc(M(:,:,zz)); axis image xy; colorbar;
    %xlabel('x'); ylabel('y'); title('Magnetization (M)');
    %%max(M(:))







return







    load charge.dat;
    load potential.dat;
    load Mx.dat;
    load My.dat;
    load Mz.dat;
    %load Hx.dat;
    %load Hy.dat;
    %load Hz.dat;

    zdim = 3;
    [ydim, xdim] = size(Mx);

    charge = reshape(charge', xdim,xdim,zdim);
    potential = reshape(potential', xdim,xdim,zdim);
    Mx = reshape(Mx', xdim,xdim,zdim);
    My = reshape(My', xdim,xdim,zdim);
    Mz = reshape(Mz', xdim,xdim,zdim);
    for z = 1:zdim
        charge(:,:,z) = charge(:,:,z)';
        potential(:,:,z) = potential(:,:,z)';
        Mx(:,:,z) = Mx(:,:,z)';
        My(:,:,z) = My(:,:,z)';
        Mz(:,:,z) = Mz(:,:,z)';
    end
    M = sqrt(Mx.^2 + My.^2 + Mz.^2);

    [ydim, xdim, zdim] = size(Mz);
    x = 0:xdim-1;
    y = 0:ydim-1;
    z = 0:zdim-1;
    [X,Y,Z] = meshgrid(x,y,z);

subplot(211);
    quiver3(X,Y,Z, Mx,My,Mz, .5);
    axis tight equal;
    xlabel('x'); ylabel('y'); title('Magnetization (M)');
    [a,b] = view();
    view(a+90-30,b);
    %view(0,90);

subplot(212);
    quiver3(X,Y,Z, Mx,My,Mz, .5);
    axis tight equal;
    xlabel('x'); ylabel('y'); title('Magnetization (M)');
    view(90,0);

    %zz = 2;
    %imagesc(M(:,:,zz)); axis image xy; colorbar;
    %xlabel('x'); ylabel('y'); title('Magnetization (M)');
    %%max(M(:))

return




%figure;
%subplot(121);
    M = sqrt(Mx.^2 + My.^2 + Mz.^2);
    imagesc(x, y, M); axis image xy;
    %imagesc(x, y, log10(M)); axis image xy;
    %caxis(log10([minM, maxM]));
    hold on;
    sf = 1;
    qh = quiver(x(1:sf:end), y(1:sf:end), Mx(1:sf:end,1:sf:end)./M(1:sf:end,1:sf:end), My(1:sf:end,1:sf:end)./M(1:sf:end,1:sf:end), 0.5);
    set(qh, 'color', 'w', 'linewidth', 1.5);
    hold off;
    xlabel('x'); ylabel('y'); title('Magnetic field (M) from FMM algorithm');
    colorbar;
%subplot(122);
    %M = sqrt(Mx.^2 + My.^2 + Mz.^2);
    %imagesc(x, y, M); axis image xy;
    %%imagesc(x, y, log10(M)); axis image xy;
    %%caxis(log10([minM, maxM]));
    %hold on;
    %sh = streamslice(x,y, Mx,My);
    %set(sh, 'color', 'w', 'linewidth', 2);
    %hold off;
    %xlabel('x'); ylabel('y'); title('Magnetic field (M) from FMM algorithm');
    %colorbar;



    %maxQ = max([charge0(:);  charge1(:); charge2(:)]);
    %minQ = min([charge0(:);  charge1(:); charge2(:)]);
    %maxV = max([potential0(:);  potential1(:); potential2(:)]);
    %minV = min([potential0(:);  potential1(:); potential2(:)]);

%figure;
%subplot(231);
    %imagesc(charge0); axis image xy; caxis([minQ, maxQ]);
%subplot(232);
    %imagesc(charge1); axis image xy; caxis([minQ, maxQ]);
%subplot(233);
    %imagesc(charge2); axis image xy; caxis([minQ, maxQ]);
%subplot(234);
    %imagesc(potential0); axis image xy; caxis([minV, maxV]);
%subplot(235);
    %imagesc(potential1); axis image xy; caxis([minV, maxV]);
%subplot(236);
    %imagesc(potential2); axis image xy; caxis([minV, maxV]);


%figure;
%subplot(121);
    %H = sqrt(Hx.^2 + Hy.^2 + Hz.^2);
    %imagesc(x, y, H); axis image xy;
    %%imagesc(x, y, log10(H)); axis image xy;
    %%caxis(log10([minH, maxH]));
    %hold on;
    %sf = 3;
    %qh = quiver(x(1:sf:end), y(1:sf:end), Hx(1:sf:end,1:sf:end)./H(1:sf:end,1:sf:end), Hy(1:sf:end,1:sf:end)./H(1:sf:end,1:sf:end), 0.5);
    %set(qh, 'color', 'w', 'linewidth', 1.5);
    %hold off;
    %xlabel('x'); ylabel('y'); title('Magnetic field (H) from FMM algorithm');
    %colorbar;
%subplot(122);
    %H = sqrt(Hx.^2 + Hy.^2 + Hz.^2);
    %imagesc(x, y, H); axis image xy;
    %%imagesc(x, y, log10(H)); axis image xy;
    %%caxis(log10([minH, maxH]));
    %hold on;
    %sh = streamslice(x,y, Hx,Hy);
    %set(sh, 'color', 'w', 'linewidth', 2);
    %hold off;
    %xlabel('x'); ylabel('y'); title('Magnetic field (H) from FMM algorithm');
    %colorbar;



return
