L1_data = load_nii('C:/ETH/Neuro/GlobalTracking/1159_L1.nii');
L2_data = load_nii('C:/ETH/Neuro/GlobalTracking/1159_L2.nii');
L3_data = load_nii('C:/ETH/Neuro/GlobalTracking/1159_L3.nii');
fa_data = load_nii('C:/ETH/Neuro/GlobalTracking/1159_FA.nii');

Lcoo = zeros(9^3,3);%numel(fa_data.img));
m=1;
for i=48-4:48+4
    for j=64-4:64+4
        for k=64-4:64+4;
Lcoo(m,1)=(L2_data.img(i,j,k)-L3_data.img(i,j,k))/(L1_data.img(i,j,k)-L3_data.img(i,j,k));
Lcoo(m,2)=L1_data.img(i,j,k)/L3_data.img(i,j,k);
            %Lcoo(m,1)=(L1_data.img(i,j,k)^2-L3_data.img(i,j,k)^2) / (L1_data.img(i,j,k)^2+L2_data.img(i,j,k)^2+L3_data.img(i,j,k)^2);
            %Lcoo(m,2)=(L1_data.img(i,j,k)^2-L2_data.img(i,j,k)^2) / (L1_data.img(i,j,k)^2+L2_data.img(i,j,k)^2+L3_data.img(i,j,k)^2);

            m=m+1;
        end
    end
end
%plot(Lcoo(:,1),0,'rx')
plot(Lcoo(:,1),Lcoo(:,2),'rx')
%plot3(Lcoo(:,1),Lcoo(:,2),Lcoo(:,3),'rx')
axis([0 1 0 10])

 snum = sum(Lcoo(:,2)<1.5)*3 + sum((Lcoo(:,2)>1.5)' .* (Lcoo(:,1)>0.5)')*2 + sum((Lcoo(:,2)>1.5)' .* (Lcoo(:,1)<0.5)')*1;

