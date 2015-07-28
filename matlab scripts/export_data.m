
nii = load_nii('C:/ETH/Neuro/Global Tractography/data_moco.nii');
dims = size(nii.img);

i=50;
j=67;
k=43;

M = zeros(1,dims(4));

            for l = 1:dims(4)
                M(l) = nii.img(i,j,k,l);
            end

csvwrite('data.csv',M);