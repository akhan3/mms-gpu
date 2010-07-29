% initial loading and setup
    clear
    meshwidth = 1e-9;
    meshdepth = 1e-9; % depth of material in z-dimension
    charge = load('charge.dat');
    potential_fmm = load('potential.dat');
    Hx_fmm = load('Hx.dat');
    Hy_fmm = load('Hy.dat');
    H_fmm = sqrt(Hx_fmm.^2 + Hy_fmm.^2);
    [ylen, xlen] = size(charge);
    x = 0:xlen-1;
    y = 0:ylen-1;
    [X, Y] = meshgrid(x,y);
    h = xlen;
    N = h^2;


% =========================================
%     Exact potential by Matlab
% =========================================

% exact potential matrix
    potential_mat = zeros(h,h);
    [yy,xx] = find(charge ~= 0);
tic;
    for k = 1:length(yy)
        fprintf('Processing %d of %d charges... %.2f%%\n', k, length(yy), k/length(yy)*100);
        potential_tmp = zeros(h,h);
        q = charge(yy(k), xx(k));
        X_ = X - (xx(k)-1);
        Y_ = Y - (yy(k)-1);
        R_ = sqrt(X_.^2 + Y_.^2);
        potential_tmp = q * 1./R_;
        potential_tmp(yy(k), xx(k)) = 0;
        %potential_mat(find(abs(potential_mat) == inf)) = 0;
        potential_mat = potential_mat + potential_tmp;
    end
time_taken = toc;

[Hx_mat, Hy_mat] = gradient(potential_mat);
constant_multiple = (meshdepth / meshwidth) / (4 * pi);
Hx_mat = constant_multiple * Hx_mat;
Hy_mat = constant_multiple * Hy_mat;
H_mat = sqrt(Hx_mat.^2 + Hy_mat.^2);


% error calculation
    abserr = abs(potential_fmm - potential_mat);
    %relerr = abserr ./ abs((potential_fmm + potential_mat)/2);
    relerr = abserr ./ abs(potential_fmm);
    relerr(isnan(relerr)) = 0;
    %relerr(find(relerr == inf)) = 0;
    max_error = max(relerr(:));
    rms_eror = sqrt(mean(relerr(:).^2));
    fprintf('Error (max = %.1e, rms = %.1e). Time taken = %g seconds\n', max_error, rms_eror, time_taken);


% find out the upper and lower limits
    maxV = max([potential_fmm(:);  potential_mat(:)]);
    minV = min([potential_fmm(:);  potential_mat(:)]);

    fh = figure;
    set(fh, 'OuterPosition', [0 0 1280 800]);
    subplot(221);
        imagesc(x,y,potential_mat); axis image xy;
        %caxis([minV, maxV]);
        colorbar;
        title('Potential from exact calculation');
    subplot(222);
        imagesc(x,y,potential_fmm); axis image xy;
        %caxis([minV, maxV]);
        colorbar;
        title('Potential from FMM algorithm');
    subplot(223);
        imagesc(x,y,charge); axis image xy;
        %spy(flipud(charge));
        colorbar;
        title('Charge distribution');
    subplot(224);
        imagesc(x,y,relerr); axis image xy; colorbar
        title_string = sprintf('Relative Error (RMS = %.1e)', rms_eror);
        title(title_string);
        %imagesc(x,y,log10(relerr)); axis image xy; colorbar;
        %title('log_{10} of Relative Error');


% find out the upper and lower limits
    maxH = max([H_fmm(:);  H_mat(:)]);
    minH = min([H_fmm(:);  H_mat(:)]);

figure;
subplot(121);
    Hx = Hx_mat;
    Hy = Hy_mat;
    H = H_mat;
    imagesc(x, y, H); axis image xy;
    hold on;
    sh = streamslice(x,y, Hx,Hy);
    set(sh, 'color', 'w');
    hold off;
    xlabel('x'); ylabel('y'); title('Magnetic field (H) from exact calculation');
    caxis([minH, maxH]);
    colorbar;
subplot(122);
    Hx = Hx_fmm;
    Hy = Hy_fmm;
    H = H_fmm;
    imagesc(x, y, H); axis image xy;
    hold on;
    sh = streamslice(x,y, Hx,Hy);
    set(sh, 'color', 'w');
    hold off;
    xlabel('x'); ylabel('y'); title('Magnetic field (H) from FMM algorithm');
    caxis([minH, maxH]);
    colorbar;

return

% compare with ohf
% =================================================
fname = ['ohf/H_',num2str(h),'x',num2str(h),'.mat'];
load(fname);
Hx_omf = Hx;
Hy_omf = Hy;
H_omf = H;
clear Hx Hy H;

% find out the upper and lower limits
    maxH = max([H_fmm(:);  H_omf(:)]);
    minH = min([H_fmm(:);  H_omf(:)]);


% error calculation
    abserr = abs(H_fmm - H_omf);
    relerr = abserr ./ abs(H_fmm);
    relerr(isnan(relerr)) = 0;
    %relerr(find(relerr == inf)) = 0;
    max_error = max(relerr(:));
    rms_eror = sqrt(mean(relerr(:).^2));
    fprintf('Error (max = %.1e, rms = %.1e).\n', max_error, rms_eror);

figure;
    %imagesc(x,y,relerr); axis image xy; colorbar
    %title_string = sprintf('Relative Error (RMS = %.1e)', rms_eror);
    imagesc(x,y,log10(relerr)); axis image xy; colorbar;
    title_string = sprintf('log_{10} of Relative Error (RMS = %.1e)', rms_eror);
    title(title_string);

figure;
subplot(121);
    Hx = Hx_omf;
    Hy = Hy_omf;
    H = H_omf;
    imagesc(x, y, log10(H)); axis image xy;
    hold on;
    sh = streamslice(x,y, Hx,Hy);
    set(sh, 'color', 'w');
    hold off;
    xlabel('x'); ylabel('y'); title('log_{10} of Magnetic field (H) from OOMMF');
    caxis(log10([minH, maxH]));
    colorbar;
subplot(122);
    Hx = Hx_fmm;
    Hy = Hy_fmm;
    H = H_fmm;
    imagesc(x, y, log10(H)); axis image xy;
    hold on;
    sh = streamslice(x,y, Hx,Hy);
    set(sh, 'color', 'w');
    hold off;
    xlabel('x'); ylabel('y'); title('log_{10} of Magnetic field (H) from FMM algorithm');
    caxis(log10([minH, maxH]));
    colorbar;



