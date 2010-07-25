% initial loading and setup
    clear
    charge = load('charge.dat');
    potential_fmm = load('potential.dat');
    [ylen, xlen] = size(charge);
    x = 0:xlen-1;
    y = 0:ylen-1;
    [X, Y] = meshgrid(x,y);
    H = xlen;
    N = H^2;

% find out the upper and lower limits
    maxV = max(potential_fmm(:));
    minV = min(potential_fmm(:));

% FMM stages animation
    %figure;
    %imagesc(x,y,charge); axis image xy; colorbar;
    %spy(flipud(charge));
    %title('Charge distribution');

    %for l = 2:log(N)/log(4)
        %filename = ['potential_L', num2str(l), '.dat'];
        %potential_fmm_L = load(filename);
        %figure;
        %imagesc(x,y,potential_fmm_L);
        %%caxis([minV, maxV]);
        %axis image xy;
        %title(['FMM at Level ', num2str(l)]);
    %end


% =========================================
%     Exact potential by Matlab
% =========================================

% exact potential matrix
    potential_exact = zeros(H,H);
    lengthBox = 10;
    meshwidth = lengthBox / sqrt(N);
    % charge setup
    [yy,xx] = find(charge ~= 0);
    for k = 1:length(yy)
        potential_tmp = zeros(H,H);
        q = charge(yy(k), xx(k));
        X_ = X - (xx(k)-1);
        Y_ = Y - (yy(k)-1);
        R_ = sqrt(X_.^2 + Y_.^2);
        potential_tmp = (q * 1./R_ / meshwidth);
        potential_tmp(yy(k), xx(k)) = 0;
        %potential_exact(find(abs(potential_exact) == inf)) = 0;
        potential_exact = potential_exact + potential_tmp;
    end

% error calculation
    relerr = abs(potential_fmm - potential_exact) ./ abs((potential_fmm + potential_exact)/2);
    relerr(isnan(relerr)) = 0;
    %relerr(find(relerr == inf)) = 0;
    max_error = max(relerr(:));
    rms_eror = sqrt(mean(relerr(:).^2));
    fprintf('Error (max = %.1e, rms = %.1e)\n', max_error, rms_eror);


% find out the upper and lower limits
    maxV = max([potential_fmm(:);  potential_exact(:)]);
    minV = min([potential_fmm(:);  potential_exact(:)]);

% plotting
    %figure;
        %imagesc(x,y,potential_exact);
        %%caxis([minV, maxV]);
        %axis image xy;
        %title('Potential from exact calculation');
    %figure;
        %%imagesc(x,y,relerr); axis image xy; colorbar
        %%title('Relative Error');
        %imagesc(x,y,log10(relerr)); axis image xy; colorbar;
        %title('log_{10} of Relative Error');

    fh = figure;
    set(fh, 'OuterPosition', [0 0 1280 800]);
    subplot(221);
        imagesc(x,y,potential_exact); axis image xy;
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



%figure;
%ph = semilogy(P, rms_error, 's-'); title(['FMM relative error (RMS)']);
%set(ph, 'linewidth',3);
%grid on
%set(gca,'XLim',[P(1)-1/2, P(end)+1/2]);
%set(gca,'XTick',P);
%set(gca,'XTickLabel',{'Monopole';'Dipole';'Quadrupole';'Octupole';'Sexdecuple'});
