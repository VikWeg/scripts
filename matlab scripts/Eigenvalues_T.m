tensor_data = load_nii('C:/ETH/Neuro/GlobalTracking/1159T/1159T_tensor.nii');
L1_data = load_nii('C:/ETH/Neuro/GlobalTracking/1159T/1159T_L1.nii');
L2_data = load_nii('C:/ETH/Neuro/GlobalTracking/1159T/1159T_L2.nii');
L3_data = load_nii('C:/ETH/Neuro/GlobalTracking/1159T/1159T_L3.nii');
fa_data = load_nii('C:/ETH/Neuro/GlobalTracking/1159T/1159T_FA.nii');
wmask = load_nii('C:/ETH/Neuro/GlobalTracking/1159T/wm_mask.nii');
wmask_diff = load_nii('C:/ETH/Neuro/GlobalTracking/1159T/wm_mask_diff.nii');

dim = size(fa_data.img);

Lcoo = zeros(sum(sum(sum(wmask_diff.img > 0))),3);%numel(fa_data.img));
m=1;
for i=1:dim(1)
    for j=1:dim(2)
        for k=1:dim(3)
            if wmask_diff.img(i,j,k) > 0
                Lcoo(m,1)=(L2_data.img(i,j,k)-L3_data.img(i,j,k))/(L1_data.img(i,j,k)-L3_data.img(i,j,k));
                Lcoo(m,2)=L1_data.img(i,j,k)/L3_data.img(i,j,k);
                m=m+1;
            end
        end
    end
end

plot(Lcoo(:,1),Lcoo(:,2),'rx')
axis([0 1 0 10])

snum = sum(Lcoo(:,2)<1.5)*3 + sum((Lcoo(:,2)>1.5)' .* (Lcoo(:,1)>0.5)')*2 + sum((Lcoo(:,2)>1.5)' .* (Lcoo(:,1)<0.5)')*1;

