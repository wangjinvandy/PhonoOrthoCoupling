clear;
%files for ROI definition

%related ROIs in orthographic visual task
spmTfiles = {};

%contrast files for estimates extraction from auditory task
ContrastFiles = {};

%volume of search
% %left hemisphere volume of search
VOI = '';

%how many subjects, praparing empty list for results
n = size(spmTfiles,1);
disp(n)
results = cell(n, 2); %preparing an array for the results
%change the directory to the one where you want to save individual ROI nii
%files
cd('')

for i=1:n
    subPath = spmTfiles{i}; %path to subjects reading spmT image
    disp(subPath)
    start = strfind(subPath, '');
    stop = strfind(subPath, '');
    subName = subPath(start:stop-2); %getting subjects code
    disp(subName)
    results{i,1} = subName;
    num = 100
    %%part of the code that creates the individual ROIs nii files,
    %%this can be run only once
    VT = spm_vol(spmTfiles{i});
    VM = spm_vol(VOI);
    
    %     T = spm_read_vols(VT).* spm_read_vols(VM);
    %     [~,I] = sort(T(:),'descend'); %This is problematic when the
    %     top-activated activaiton is negative corrected by Jin 2/21/2022
    %     R = zeros(size(T));
    %     R(I(1:100)) = 1;
    
    %%%%%%%%%%%%%added by Jin%%%%%%%%%%%%%%%%%%%%
    m=spm_read_vols(spm_vol(VM));
    mask_indx=find(m>0);
    [Tdata,~]=spm_read_vols(VT);
    index_non0=find(Tdata(mask_indx)~=0); %this step cleans the 0-activated voxels in the vOT mask
    non0_mask_indx=mask_indx(index_non0); 
    values_inmask=Tdata(non0_mask_indx);
    [top_activation,I]=sort(values_inmask,'descend'); %this step only sort the non-zero values in the mask
    top_indx=non0_mask_indx(I(1:num)); 
    R = zeros(size(Tdata));
    R(top_indx) = 1;
    %%%%%%%%%%%%%%end of addition%%%%%%%%%%%%%%%%%
    
    VM.fname = strcat(subName,'.nii');
    spm_write_vol(VM,R);
    %%end of the creating ROI file part
    ConEst = spm_summarise(ContrastFiles{i},strcat(subName,'.nii'),'', true);
    results{i,2} = mean(ConEst, 'omitnan');
end
% Table = cell2table(results, 'VariableNames', {'Code', 'MeanBeta'});
% writetable(Table, '....csv', 'Delimiter', ',');
