function [region,percent] = find_roi(vol_coordinate,aparc_vol,aparc_aseg,aparc_indx,radius)
%
%
%
%

if ndims(aparc_vol) == 3
    %% roi calculation for 3D nifti template
    for icoors = 1:size(vol_coordinate,1)
        current_coor = round(vol_coordinate(icoors,:));
        n = 1;
        all_region_indx = [];
        
        for xt = current_coor(1)-radius:current_coor(1)+radius
            for yt = current_coor(2)-radius:current_coor(2)+radius
                for zt = current_coor(3)-radius:current_coor(3)+radius
                    if norm(current_coor-[xt,yt,zt]) < radius
                        all_region_indx(n) = aparc_vol(xt,yt,zt);
                        n = n+1;
                    end
                end
            end
        end
        if ~any(all_region_indx)
            region{icoors} = 'NO_ROI';
            percent{icoors} = 0;
        else
            if sum(all_region_indx>0) > round(length(all_region_indx)/2)
                all_region_indx_clear = all_region_indx;
                all_region_indx_clear(all_region_indx_clear==0) = [];
                [mode_row, occu_num] = mode(all_region_indx_clear');
                percent{icoors} = occu_num/size(all_region_indx,2);
                region{icoors} = deblank(aparc_aseg(aparc_indx == mode_row,:));
            else
                region{icoors} = 'NO_ROI';
                percent{icoors} = 0;
            end
        end
        
    end
    
    
else
    %% roi calculation for 4D nifti template
    % resize coordinates to fit template size
    vol_size = size(aparc_vol);
    resize_mat = eye(3);
    resize_mat(1) = vol_size(1)/256;
    resize_mat(5) = vol_size(2)/256;
    resize_mat(9) = vol_size(3)/256;
    vol_coordinate = vol_coordinate*resize_mat;
    
    if size(aparc_aseg,1) == size(aparc_indx,1)
        
        for icoors = 1:size(vol_coordinate,1)
            current_coor = round(vol_coordinate(icoors,:));
            max_roi_indx = [];
            max_prob = 0;
            for iroi = 1:size(aparc_indx,1)
                aparc_vol_single = aparc_vol(:,:,:,iroi);
                if aparc_vol_single(current_coor(1),current_coor(2),current_coor(3)) > 0 &&  ...
                        aparc_vol_single(current_coor(1),current_coor(2),current_coor(3)) > max_prob;
                    max_prob = aparc_vol_single(current_coor(1),current_coor(2),current_coor(3));
                    max_roi_indx  = iroi;
                end
            end
            
            if isempty(max_roi_indx)
                region{icoors} = 'NO_ROI';                
                percent{icoors} = 0;
            else
                percent{icoors} = max_prob/100;
                region{icoors} = deblank(aparc_aseg(max_roi_indx,:));
            end
            
        end
        
    else 
        errordlg('Unmatched template labeling file!\n')
    end
    

    
end
percent = percent';
region = region';
    
end


