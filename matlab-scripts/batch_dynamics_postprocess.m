function batch_dynamics_postprocess(folder_name)


numfiles = 50;
Istart = 1e-6;  %Starting current in A
Iend = 5e-3;    %End current in A

Ivec = linspace(Istart,Iend,numfiles);
Nmax = 10006;
% Ivec = logspace(log10(Istart), log10(Iend), numfiles);
% Nmax = 19994;

% declare arrays
time = zeros(Nmax,50);
tindex = zeros(size(time));
Mx = zeros(size(time));

% loop for loading files
for fn = 1:numfiles
    dirName = [folder_name,'/sim_STO', num2str(fn,'%.4d')];
%     dirName = [folder_name,'/STO', num2str(fn,'%.4d')];
    fileName = [dirName,'/dynamics.dat'];
    disp(fileName)
    dynamics = load(fileName);
    n = length(dynamics);

    % fill the arrays
    tindex(1:n,fn) = dynamics(:,1);
    time(1:n,fn) = dynamics(:,2);
    Mx(1:n,fn) = dynamics(:,5);
    % stretch the array
    if n ~= Nmax
        Mx(n+1:end,fn) = Mx(n,fn) * ones(size(Mx(n+1:end,fn)));
        dt = time(2,fn) - time(1,fn);
        time(n+1:end,fn) = time(n,fn) + (1:Nmax-n) * dt;
    end

end

save dynamics_pp

N1 = 1;
N2 = 50;
plot(time(:,N1:N2)/1e-9, Mx(:,N1:N2))

