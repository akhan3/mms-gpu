% initial loading and setup
    clear

    q = 1.6e-19;
    load dynamics.dat
    index = dynamics(:,1);
    t = dynamics(:,2);
    energy = dynamics(:,3);
    Mx = dynamics(:,4);
    My = dynamics(:,5);
    Mz = dynamics(:,6);
    M = dynamics(:,7);

%set(gcf, 'OuterPosition', [0 0 1280 800]);
subplot(221);
    plot(t, energy/q, '.');
    xlabel('time'); ylabel('Energy (eV)');

subplot(222);
    plot(t, Mx, t, My, t, Mz, t, M);
    legend('Mx', 'My', 'Mz', 'M');
    xlabel('time'); ylabel('Magnetization (A/m)');


% load and post process M files
    load Mx.dat;
    load My.dat;
    load Mz.dat;
    [ydim, xdim] = size(Mx);
    zdim = 1;
    Mx = reshape(Mx', xdim,xdim,zdim);
    My = reshape(My', xdim,xdim,zdim);
    Mz = reshape(Mz', xdim,xdim,zdim);
    for z = 1:zdim
        Mx(:,:,z) = Mx(:,:,z)';
        My(:,:,z) = My(:,:,z)';
        Mz(:,:,z) = Mz(:,:,z)';
    end
    M = sqrt(Mx.^2 + My.^2 + Mz.^2);
    [ydim, xdim] = size(Mx);
    [ydim, xdim, zdim] = size(Mz);
    x = 0:xdim-1;
    y = 0:ydim-1;
    z = 0:zdim-1;
    [X,Y,Z] = meshgrid(x,y,z);


%subplot(223);
    %quiver3(X,Y,Z, Mx,My,Mz, .5);
    %axis tight equal;
    %xlabel('x'); ylabel('y'); title('Magnetization (M)');
    %[a,b] = view();
    %view(a+90-30,b);
    %%view(0,90);

%subplot(224);
    %quiver3(X,Y,Z, Mx,My,Mz, .5);
    %axis tight equal;
    %xlabel('x'); ylabel('y'); title('Magnetization (M)');
    %view(90,0);



% load and post process H files
    load Hx.dat;
    load Hy.dat;
    load Hz.dat;
    [ydim, xdim] = size(Hx);
    %zdim = xdim;
    Hx = reshape(Hx', xdim,xdim,zdim);
    Hy = reshape(Hy', xdim,xdim,zdim);
    Hz = reshape(Hz', xdim,xdim,zdim);
    for z = 1:zdim
        Hx(:,:,z) = Hx(:,:,z)';
        Hy(:,:,z) = Hy(:,:,z)';
        Hz(:,:,z) = Hz(:,:,z)';
    end
    H = sqrt(Hx.^2 + Hy.^2 + Hz.^2);
    [ydim, xdim] = size(Hx);
    [ydim, xdim, zdim] = size(Hz);
    x = 0:xdim-1;
    y = 0:ydim-1;
    z = 0:zdim-1;
    [X,Y,Z] = meshgrid(x,y,z);

return

%load fmm_data_64x64.mat
% draw H field in slices
    %figure('Colormap',map)
    slice_num = zdim/2;
    %contourslice(H,[],[],slice_num)
    contourslice(H,[],[],[1, zdim/2, zdim])
    axis image xy
    view(3)
    xlabel('x'); ylabel('y'); zlabel('z');


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
