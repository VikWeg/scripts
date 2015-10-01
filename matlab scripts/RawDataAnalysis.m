L1 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/L1.nii');
L2 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/L2.nii');
L3 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/L3.nii');
T = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/tensor.nii');
wmask = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/wm_mask.nii');

dimL = size(L1.img);
sizeL = sum(sum(sum(wmask.img > 0)));
dimT = size(T.img);
sizeT = numel(T.img);

DX = (L2.img - L3.img)./(L1.img - L3.img);
DXm = DX( wmask.img > 0 );
DiffPlaneX = reshape(DXm,[1,numel(DXm)]);

DY = (L1.img)./(L3.img);
DYm = DY( wmask.img > 0 );
DiffPlaneY = reshape(DYm,[1,numel(DYm)]);

Valid = find(  (DiffPlaneY < 10) .* (DiffPlaneY > 1) .* (DiffPlaneX > 0) .* (DiffPlaneX < 1) );

DiffPlaneX = DiffPlaneX(Valid);
DiffPlaneY = DiffPlaneY(Valid);

% subplot(1,2,1)
% plot(DiffPlaneX,DiffPlaneY,'rx')
% axis square

DiffData = [DiffPlaneX' -DiffPlaneY'];

DiffHist = hist3(DiffData,[20 20])/length(DiffData)*100;

% subplot(1,2,2)
imagesc([0 1],[1 10],DiffHist')
colorbar
axis square
