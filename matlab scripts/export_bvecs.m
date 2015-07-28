bvecsfile=fopen('bvecs.csv','w+');
bvecs = importdata('C:/ETH/Neuro/Global Tractography/ep2d_diff_rolled.bvecs');


for i = 1:length(bvecs)
    fprintf(bvecsfile,'%d,%d,%d\n',bvecs(i,1:3));  
end

fclose(bvecsfile);