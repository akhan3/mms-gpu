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

%% FMM stages animation
    %figure; 
    %imagesc(x,y,charge); axis image xy; colorbar;
    %%spy(flipud(charge));
    %title('Charge distribution');
    %for l = 2:log(N)/log(4)
    %%for l = log(N)/log(4):log(N)/log(4)
        %filename = ['potential_L', num2str(l), '.dat'];
        %potential_fmm_L = load(filename); 
        %figure; 
        %imagesc(x,y,potential_fmm_L); 
        %caxis([minV, maxV]);
        %axis image xy;
        %title(['FMM at Level ', num2str(l)]);
    %end


% =========================================
%     Exact potential by Matlab
% =========================================

% exact potential matrix
    potential_exact = zeros(H,H);
    % charge setup
    [yy,xx] = find(charge ~= 0);
    for k = 1:length(yy)
        q = charge(yy(k), xx(k));
        X_ = X - (xx(k)-1);
        Y_ = Y - (yy(k)-1);
        R_ = sqrt(X_.^2 + Y_.^2);
        potential_exact = potential_exact + (q * 1./R_);
        potential_exact(find(abs(potential_exact) == inf)) = 0;
    end

% error calculation
    relerr = abs(potential_fmm - potential_exact) ./ abs(potential_exact);
    relerr(isnan(relerr)) = 0;
    relerr(find(relerr == inf)) = 0;
    max_error = max(relerr(:))
    rms_eror = sqrt(mean(relerr(:).^2))
    

% find out the upper and lower limits
    maxV = max([potential_fmm(:);  potential_exact(:)]);
    minV = min([potential_fmm(:);  potential_exact(:)]);

% plotting
    figure;
    subplot(221);
        %imagesc(x,y,charge); axis image xy;
        spy(flipud(charge));
        title('Charge distribution');
    subplot(222);
        imagesc(x,y,potential_fmm); axis image xy;
        caxis([minV, maxV]);
        colorbar;
        title('Potential from FMM algorithm');
    subplot(223);
        imagesc(x,y,potential_exact); axis image xy;
        caxis([minV, maxV]);
        colorbar;
        title('Potential from exact calculation');
    subplot(224);
        %imagesc(x,y,relerr); axis image xy; colorbar
        imagesc(x,y,log10(relerr)); axis image xy; colorbar;
        title('log_{10} of relative error');
