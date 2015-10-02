%% Load
L1 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/L1.nii');
L2 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/L2.nii');
L3 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/L3.nii');
T = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/tensor.nii');
wmask = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/wm_mask.nii');

%% Preselection
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

Valid = find(  (DiffPlaneY < 7) .* (DiffPlaneY > 1) .* (DiffPlaneX > 0) .* (DiffPlaneX < 1) );

%% Data
DiffPlaneX = DiffPlaneX(Valid);
DiffPlaneY = DiffPlaneY(Valid);

DiffData = [DiffPlaneX' DiffPlaneY'];
DiffHist = hist3(DiffData,[20 20])/length(DiffData)*100;

%% Plot

lim1 = 1;
lim2 = 2;
steps = 20;

sphere = zeros(1,steps);
cigar = zeros(1,steps);
discus = zeros(1,steps);
tau = zeros(1,steps);
count = length(DiffPlaneY);

for t = 1:steps
    tau(t) = lim1 + (lim2-lim1)/(steps-1)*(t-1);
    sphere(t)= sum(  DiffPlaneY < tau(t))/count;
    cigar(t) = sum( (DiffPlaneY > tau(t)) .* (DiffPlaneX < 0.5) )/count;
    discus(t)= sum( (DiffPlaneY > tau(t)) .* (DiffPlaneX > 0.5) )/count;
end

% plot(tau,sphere,'ro',tau,cigar,'bd',tau,discus,'gh')
% legend('sphere','cigar','discus')

subplot(1,2,1)
plot(DiffPlaneX,DiffPlaneY,'')
axis([0 1 1 7])
hold on
imagesc([0 1],[1 7],DiffHist')
colorbar
hold off

subplot(1,2,2)
plot(DiffPlaneX,DiffPlaneY,'w')
axis([1 lim2 0 1])
hold on
bar(tau,[cigar' discus' sphere'],'stacked')
hold off
%%




















