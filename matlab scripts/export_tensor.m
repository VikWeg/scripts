%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subject = '1159T';
size = [40 94 58];
center = [85 67 56];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tensor = load_nii(['C:/ETH/Neuro/GlobalTracking/subjects/' subject '/tensor.nii']);

[tensorfile,message] = fopen(['C:/ETH/Neuro/GlobalTracking/subjects/' subject ...
    '/tensor_c_' num2str(center(1)) '_' num2str(center(2)) '_' num2str(center(3)) ...
    '_s_' num2str(size(1)) '_' num2str(size(2)) '_' num2str(size(3)) '.txt'],'w');

fprintf(tensorfile,'{');
for i = center(1)-floor(size(1)/2):center(1)+floor(size(1)/2)
    fprintf(tensorfile,'{');
    for j = center(2)-floor(size(2)/2):center(2)+floor(size(2)/2)
        fprintf(tensorfile,'{');
        for k = center(3)-floor(size(3)/2):center(3)+floor(size(3)/2)
            fprintf(tensorfile,'{');
            for t = 1:6
                if t < 6
                fprintf(tensorfile,'%f,',tensor.img(i,j,k,t));
                else
                fprintf(tensorfile,'%f',tensor.img(i,j,k,t));
                end
            end
            if k < center(3)+floor(size(3)/2)
            fprintf(tensorfile,'},');
            else
            fprintf(tensorfile,'}');
            end
        end
        if j < center(2)+floor(size(2)/2)
        fprintf(tensorfile,'},');
        else
        fprintf(tensorfile,'}');
        end
    end
    if i < center(1)+floor(size(1)/2)
    fprintf(tensorfile,'},');
    else
    fprintf(tensorfile,'}');
    end
end
fprintf(tensorfile,'}');

fclose(tensorfile);