

function coordinates = ccep_CreateROIs(subinfo)
% Create ROI based on the electrode coordinates
filename = 'elec_image.nii';
temInfo = load_untouch_nii([subinfo.mainpath filesep 'brain3D' filesep 'T1.nii']);
if exist([subinfo.mainpath filesep 'brain3D' filesep 'coordinates.txt' ],'file')
    coordinates = dlmread([subinfo.mainpath filesep 'brain3D' filesep 'coordinates.txt']);
elseif exist([subinfo.mainpath filesep 'brain3D' filesep 'autocoordinates.mat' ],'file')
    coordinates = load([subinfo.mainpath filesep 'brain3D' filesep 'autocoordinates.mat']);
    coordinates = coordinates.savecoors(:,end-2:end);
end

if ~exist([subinfo.mainpath filesep 'brain3D' filesep filename],'file')
   
    % How many rois do they want to drop into this mask file.
    
    nROIS = size(coordinates,1);
    
    if nROIS < 1
        errordlg('Failed to load electrode coordinates file.')
        return
    end
    
    % Get the name of the reference image.
    
    
    vox2ras = [temInfo.hdr.hist.srow_x;temInfo.hdr.hist.srow_y;temInfo.hdr.hist.srow_z;0 0 0 1];
    
    dim = size(temInfo.img);
    
    maskIMG = zeros(size(temInfo.img));
    
    
    nii.hdr = temInfo.hdr;
    
    roiINFO = {};
    
    % Now loop and get the information for the ROI.
    
    for iROI = 1:nROIS
        roiINFO{iROI}.center_mm = coordinates(iROI,:)';
        
        tmp = inv(vox2ras)*([roiINFO{iROI}.center_mm ;1]);
        roiINFO{iROI}.center_vox = tmp(1:3);
        roiINFO{iROI}.type = 'Sphere';% 'Box'
        switch roiINFO{iROI}.type
            case 'Sphere'
                str = 'Radius (mm)';
                nIn = 1;
            case 'Box'
                str = 'Straddle Size (mm)';
                nIn = 3;
        end
        roiINFO{iROI}.size = nIn;
    end
    
    % Build an array of coordinates for each and every voxel in the mask
    
    % PMat = spm_imatrix(P_HDR.mat);
    % xSize = PMat(7);
    % ySize = PMat(8);
    % zSize = PMat(9);
    
    xOrds = (1:dim(1))'*ones(1,dim(2));
    yOrds = ones(dim(1),1)*(1:dim(2));
    xOrds = xOrds(:)';
    yOrds = yOrds(:)';
    
    Coords = zeros(3,prod(dim(1:3)));
    
    for iZ = 1:dim(3)
        zOrds = iZ*ones(1,length(xOrds));
        Coords(:,(iZ-1)*length(xOrds)+1:iZ*length(xOrds)) = [xOrds; yOrds; zOrds];
    end
    
    % Now put them into mm's
    
    mmCoords = vox2ras*[Coords;ones(1,size(Coords,2))];
    mmCoords = mmCoords(1:3,:);
    
    boxBIT = zeros(4,size(mmCoords,2));
    
    % Now loop on the ROI definitions and drop them
    % into the mask image volume matrix.
    
    for iROI = 1:nROIS
        % Found the center of this ROI in voxels
        % and then build it.
        xs = mmCoords(1,:) - roiINFO{iROI}.center_mm(1);
        ys = mmCoords(2,:) - roiINFO{iROI}.center_mm(2);
        zs = mmCoords(3,:) - roiINFO{iROI}.center_mm(3);
        switch roiINFO{iROI}.type
            case 'Sphere'
                radii = sqrt(xs.^2+ys.^2+zs.^2);
                VOXIdx = find(radii<=roiINFO{iROI}.size);
            case 'Box'
                xsIDX = find(abs(xs)<=roiINFO{iROI}.size(1));
                ysIDX = find(abs(ys)<=roiINFO{iROI}.size(2));
                zsIDX = find(abs(zs)<=roiINFO{iROI}.size(3));
                boxBIT  = 0*boxBIT;
                boxBIT(1,xsIDX) = 1;
                boxBIT(2,ysIDX) = 1;
                boxBIT(3,zsIDX) = 1;
                boxBIT(4,:) = boxBIT(1,:).*boxBIT(2,:).*boxBIT(3,:);
                VOXIdx = find(boxBIT(4,:));
        end
        maskIMG(VOXIdx) = 1000;
    end
    nii.img = maskIMG;
    % Now write out the image.
    
    save_nii(nii, [subinfo.mainpath filesep 'brain3D' filesep filename])
    
    fprintf('\nFinished Building ROI Image.\n');
end

%
% All done.
%
