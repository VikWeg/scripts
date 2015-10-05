%% Load
L1 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/L1.nii');
L2 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/L2.nii');
L3 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/L3.nii');
T = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/tensor.nii');
wmask = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/wm_mask.nii');

L1_2 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/1159T_L1.nii');
L2_2 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/1159T_L2.nii');
L3_2 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/1159T_L3.nii');
T_2 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/1159T_tensor.nii');
wmask_2 = load_nii('C:/ETH/Neuro/GlobalTracking/subjects/1159T/wm_mask_diff.nii');

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

subplot(4,4,3)
plot(DiffPlaneX,DiffPlaneY,'x')

Valid = find(  (DiffPlaneY < 7) .* (DiffPlaneY > 1) .* (DiffPlaneX > 0) .* (DiffPlaneX < 1) );

%% Data
DiffPlaneX = DiffPlaneX(Valid);
DiffPlaneY = DiffPlaneY(Valid);

DiffData = [DiffPlaneX' DiffPlaneY'];
DiffHist = hist3(DiffData,[40 40])/length(DiffData)*100;

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

subplot(4,4,1)
plot(DiffPlaneX,DiffPlaneY,'')
axis([0 1 1 7])
hold on
imagesc([0 1],[1 7],DiffHist')
colorbar
hold off

subplot(4,4,2)
plot(DiffPlaneX,DiffPlaneY,'w')
axis([1 lim2 0 1])
hold on
bar(tau,[cigar' discus' sphere'],'stacked')
hold off
%% Eigen Distribution
L1vec = reshape(L1.img(:,:,:),[1,numel(L1.img)]);
L2vec = reshape(L1.img(:,:,:),[1,numel(L1.img)]);
L3vec = reshape(L3.img(:,:,:),[1,numel(L1.img)]);

subplot(4,4,4)
histogram(L1vec,linspace(0.0002,0.0017))
title('L1')

subplot(4,4,5)
histogram(L2vec,linspace(0.0002,0.0017))
title('L2')

subplot(4,4,6)
histogram(L3vec,linspace(0.0002,0.0017))
title('L3')

T1vec = reshape(T.img(:,:,:,1),[1,numel(L1.img)]);
T2vec = reshape(T.img(:,:,:,2),[1,numel(L1.img)]);
T3vec = reshape(T.img(:,:,:,3),[1,numel(L1.img)]);
T4vec = reshape(T.img(:,:,:,4),[1,numel(L1.img)]);
T5vec = reshape(T.img(:,:,:,5),[1,numel(L1.img)]);
T6vec = reshape(T.img(:,:,:,6),[1,numel(L1.img)]);

subplot(4,4,7)
histogram(T1vec,linspace(0.0002,0.0017))
title('T1')

subplot(4,4,8)
histogram(T2vec,linspace(0.00005,0.0005))
title('T2')

subplot(4,4,9)
histogram(T3vec,linspace(0.00005,0.0005))
title('T3')

subplot(4,4,10)
histogram(T4vec,linspace(0.00005,0.002))
title('T4')

subplot(4,4,11)
histogram(T5vec,linspace(0.00005,0.0005))
title('T5')

subplot(4,4,12)
histogram(T6vec,linspace(0.00005,0.002))
title('T6')

%% Subject Comparison

L1_2vec = reshape(L1_2.img(:,:,:),[1,numel(L1_2.img)]);

subplot(4,4,13)
plot(L1vec,L1_2vec,'x',[-0.001 0.002],[-0.001 0.002],'r')
title('L1 vs L1_2 abs')

relErr1 = abs(L1vec - L1_2vec)./L1vec;

sub = randperm(length(L1vec));
sub = sub(1:ceil(length(L1vec)/10));

subplot(4,4,14)
plot(L1vec(sub),relErr1(sub),'x')
axis([0.0001 0.002 0 1])
title('L1 vs L1_2 rel')

%% L1 Err Hist
val = find( (L1vec > 0.0001) .* (relErr1 < 1) );

VSHist = hist3([L1vec(val)' relErr1(val)'],[20 20])/length(val)*100;

%subplot(4,4,15)
plot(0.0001,0,'')
axis([0.0001 0.002 0 1])
hold on
imagesc([0.0001 0.002],[0 1],VSHist')
colorbar
hold off

%% L2 Err Hist

L2_2vec = reshape(L2_2.img(:,:,:),[1,numel(L2_2.img)]);

relErrL2 = abs(L2vec - L2_2vec)./L2vec;

valL2 = find( (L2vec > 0.0001) .* (relErrL2 < 1) );

ErrL2Hist = hist3([L2vec(valL2)' relErrL2(valL2)'],[20 20])/length(valL2)*100;

%subplot(4,4,16)
plot(0.0001,0,'')
axis([0.0001 0.002 0 1])
hold on
imagesc([0.0001 0.002],[0 1],ErrL2Hist')
colorbar
hold off

%% L3 Err Hist

L3_2vec = reshape(L3_2.img(:,:,:),[1,numel(L3_2.img)]);

relErrL3 = abs(L3vec - L3_2vec)./L3vec;

valL3 = find( (L3vec > 0.0001) .* (relErrL3 < 1) );

ErrL3Hist = hist3([L3vec(valL3)' relErrL3(valL3)'],[20 20])/length(valL3)*100;

%subplot(4,4,16)
plot(0.0001,0,'')
axis([0.0001 0.002 0 1])
hold on
imagesc([0.0001 0.002],[0 1],ErrL3Hist')
colorbar
hold off








