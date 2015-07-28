fa_data = load_nii('C:/ETH/Neuro/GlobalTracking/1159_FA.nii');

FA = zeros(1,9^3);%numel(fa_data.img));
m=1;
for i=64-4:64+4
    for j=64-4:64+4
        for k=48-4:48+4;
            FA(m)=fa.img(i,j,k);
            m=m+1;
        end
    end
end