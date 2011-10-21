function Minit_final(folder_name)

% figure;
[X Y Z] = meshgrid(1:13, 1:13, 1);

subplot(121);
Mi = load([folder_name, '/Minit.dat']);
Mi = reshape(Mi',3,13,13);
% quiver(1:13,1:13, squeeze(Mi(1,:,:)), squeeze(Mi(2,:,:)) ); axis image;
quiver3(X,Y,Z, squeeze(Mi(1,:,:)), squeeze(Mi(2,:,:)), squeeze(Mi(3,:,:)) ); 
axis image; view(2); grid off;


subplot(122);
Mf = load([folder_name, '/Mfinal.dat']);
Mf = reshape(Mf',3,13,13);
% quiver(1:13,1:13, squeeze(Mf(1,:,:)), squeeze(Mf(2,:,:)) ); axis image;
quiver3(X,Y,Z, squeeze(Mf(1,:,:)), squeeze(Mf(2,:,:)), squeeze(Mf(3,:,:)) ); 
axis image; view(2); grid off;
