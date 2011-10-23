%This creates a graph of Fourier transformed signals

Makeplots = 1;

if Makeplots
    figure;
    plot(time(:,1:5:end)/1e-9, Mx(:,1:5:end))
end
xlabel('t [ns]')
ylabel('m_x(t)');


figure; hold on;
for i = 1:numfiles
    Mxsignal = Mx(:,i);
    tsignal = time(:,i);

% Now comes the Fourier transformation part
    svec = Mxsignal - mean(Mxsignal);
    dt = tsignal(2)-tsignal(1);
    fsample = 1/dt;
    NFFT = 2^nextpow2(length(svec));
    Y = fft(svec,NFFT)/length(svec);
    fvec(i,:) = fsample/2*linspace(0,1,NFFT/2+1);
    Yvec(i,:) = 2*abs(Y(1:NFFT/2+1));
%     plot(fvec(i,:),Yvec(i,:))
end

plot(fvec(:,:)'/1e9, Yvec(:,:)')
title('Single-Sided Amplitude Spectrum of svec(t)')
xlabel('Frequency (GHz)')
ylabel('|Y(f)|');
set(gca, 'xlim', [0 5]);

%it is now converting the entire image into a contour plot
%Shows the Fourier transform amplitude as a function of current

fpoints = 1000;
%fmin = min(min(fvec));
%fmax = max(max(fvec));

fmin = 0.0;
fmax = 10.0E+9;


fplot = linspace(fmin,fmax,fpoints);

FFTmx(1:numfiles,1:fpoints) = 0.0;

for i = 1 : numfiles
    if (1)
        finterp = interp1(fvec(i,:),Yvec(i,:),fplot);
        FFTmx(i,:) = finterp;
    end
end

figure;
alpha = 1.;
surf(Ivec/1e-3, fplot/1e9, (transpose(FFTmx).^alpha));
shading flat
axis xy
set(gca, 'ylim', [0 5]);
xlabel('I [mA]')
ylabel('Frequency [GHz]');
title('Signal Power is denoted by the color');
colorbar
colormap jet

