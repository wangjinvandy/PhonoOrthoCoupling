%%%%correlational MVPA analysis, written by Jin 2/21/2022 %%%%

clear;
%files for ROI definition


%volume of search
VOI_path='';
VOI=''; 

% Visual spmT contrast of interest
timepoint='tp1'; %'tp1'
spmT_visual_path=''; spmT_visual=''; 
% Phono spmT contrast of interest
spmT_audio_path='\'; spmT_audio= '';

subjects= {
'sub-Dys001'

}; %your participant number (put all participants)
num=100; %number of top voxels you want to select

output_path = ''
%output file 
writefile='name';
cd(output_path);
if exist(writefile)
   delete(writefile);
end
fid_w=fopen(writefile,'w');
fprintf(fid_w, '%s %4s %.4s %4s\n', 'participant_id', 'r_correl','p_correl','rz_correl');


for i=1:length(subjects)
    %cd([root, filesep, subjects{i}]);
    if timepoint=='tp2'
    spmTfile_1 = [spmT_visual_path,filesep,subjects{i},filesep, 'output_...', filesep, spmT_visual,'.nii']; %spmT_0031 visual_task words_falsefonts_tp2
    else
    spmTfile_1 = [spmT_visual_path,filesep,subjects{i},filesep, 'output_...', filesep, spmT_visual,'.nii']; %spmT_0031 visual_task words_falsefonts_tp1
    end
    spmTfile_2 = [spmT_audio_path,filesep, subjects{i},filesep, 'output_....', filesep, spmT_audio,'.nii']; %spmT_0039 aliteration_tp2
    VOIfile=[VOI_path,filesep,VOI,'.nii'];
    VM = spm_vol(VOIfile);
    VT = spm_vol(spmTfile_1);
    
    %     T = spm_read_vols(VT).* spm_read_vols(VM);
    %     [x,I] = sort(T(:),'descend'); %this will cause a problem when there
    %     is negative activaiton in the top activated voxels, detected by Jin
    %     R = zeros(size(T));
    %     R(I(1:num)) = 1;
    
    %%added by Jin Wang, 2/21/2021%%%%%%%%%%
    m=spm_read_vols(spm_vol(VM));
    mask_indx=find(m>0);
    [Tdata,~]=spm_read_vols(VT);
    index_non0=find(Tdata(mask_indx)~=0);
    non0_mask_indx=mask_indx(index_non0);
    values_inmask=Tdata(non0_mask_indx);
    [top_activation,I]=sort(values_inmask,'descend');
    top_indx=non0_mask_indx(I(1:num));
    R = zeros(size(Tdata));
    R(top_indx) = 1;
    %%%%%%%%%%%end of modification %%%%%%%%%%%%
    
%     VM.fname = strcat(subjects{i},'_', VOI,'_top-',num2str(num), '_', spmT_visual,'_',timepoint, '.nii');
%     new_mask=spm_write_vol(VM,R); 
    %end of creating the top-activated voxels ROI
    
    %extracting the t-values from spmT maps at the top voxels' location
    [x,y,z]=ind2sub(size(R),top_indx);
    XYZ=[x,y,z]';
    value_spmT1=spm_get_data(spmTfile_1,XYZ);
    value_spmT2=spm_get_data(spmTfile_2,XYZ);
    [correl,p]=corrcoef(value_spmT1,value_spmT2); %calculate the correlations between the two spmT maps
    r_correl=correl(1,2);
    p_correl=p(1,2);
    rz_correl=0.5*log((1+r_correl)/(1-r_correl)); %conduct the fisher-z transformation on the correlational values
    
    fprintf(fid_w,'%s %4f %4f %4f',subjects{i},r_correl,p_correl,rz_correl);
    fprintf(fid_w,'\n');
    
end


