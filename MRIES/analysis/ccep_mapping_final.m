function [n,stiLocation] = ccep_mapping_final(subinfo,value,label)

warning('off')
files = dir([subinfo.mainpath filesep 'data' filesep 'ccep*.mat']);
elec1 = [];
elec2 = [];
for j =1:length(files)
    [~, filename, ext] = fileparts(files(j).name);
    indx = strfind(filename,'_');
    elec1 = [elec1;str2double(filename(indx(end-1)+1:indx(end)-1))]; %ok<ST2NM>
    elec2 = [elec2;str2double(filename(indx(end)+1:end))]; 
end
elec1 = sort(elec1); elec2 = sort(elec2);

temInfo = load_untouch_nii([subinfo.mainpath filesep 'brain3D' filesep 'T1.nii']);
if exist([subinfo.mainpath filesep 'brain3D' filesep 'coordinates.txt' ],'file')
    coordinates = dlmread([subinfo.mainpath filesep 'brain3D' filesep 'coordinates.txt']);
elseif exist([subinfo.mainpath filesep 'brain3D' filesep 'autocoordinates.mat' ],'file')
    coordinates = load([subinfo.mainpath filesep 'brain3D' filesep 'autocoordinates.mat']);
    coordinates = coordinates.savecoors(:,end-2:end);
end
    
vox2ras = [temInfo.hdr.hist.srow_x;temInfo.hdr.hist.srow_y;temInfo.hdr.hist.srow_z;0 0 0 1];
elecpos = vox2ras\[coordinates ones(size(coordinates,1),1)]';
elecpos = round(elecpos(1:3,:)');
stiLocation=elecpos([label(1),label(2)],:);
data = [elecpos value];
ndim = size(temInfo.img);

method = 'linear+nearest';
if strcmp(method,'gauss')
    
    [data_interpolated] = ccep_interpolate_gauss(ndim,data);
    n = round( temInfo.hdr.hist.srow_x(1));
elseif strcmp(method,'linear+nearest')
    
    n = 2;
    %     bb = [ round(temInfo.hdr.hist.qoffset_x)  round(temInfo.hdr.hist.qoffset_y) round(temInfo.hdr.hist.qoffset_z)
    %         temInfo.hdr.dime.dim(2)-abs(round(temInfo.hdr.hist.qoffset_x))-1   temInfo.hdr.dime.dim(3)-abs(round(temInfo.hdr.hist.qoffset_y))-1 temInfo.hdr.dime.dim(4)-abs(round(temInfo.hdr.hist.qoffset_z))-1 ];
    % interpolate
    dmax1 = 5;  % changed by Liang 2017/04/23
    [data_interpolated,index1] = ccep_interpolate_ln(ndim,data,dmax1,n,'nearest');
    dmax2 = 10;  % changed by Liang 2017/04/23
    [di2,index2] = ccep_interpolate_ln(ndim,data,dmax2,n,'linear');
    index3 = intersect(find(isnan(data_interpolated)),index2);
    data_interpolated(index3) = di2(index3);
%     data_interpolated=permute(data_interpolated,[2 1 3]);
    temInfo.hdr.dime.dim(2:4) = size(data_interpolated);
end


temInfo.img = data_interpolated;
temInfo.hdr.dime.datatype = 4;
temInfo.hdr.dime.glmax = max(data_interpolated(:));
temInfo.hdr.dime.glmin = min(data_interpolated(:));
save_untouch_nii(temInfo,[subinfo.mainpath filesep 'brain3D' filesep 'ccep_map_raw.nii'])

matlabbatch{1}.spm.spatial.smooth.data = {[[subinfo.mainpath filesep 'brain3D' filesep 'ccep_map_raw.nii'] ',1']};
matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 1;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm('defaults', 'EEG');
spm_jobman('run', matlabbatch);


end


function [data_interpolated] = ccep_interpolate_gauss(ndim,data)

data = single(data);
data_interpolated = zeros(ndim,'single');
aind = 1:prod(ndim);
[i,j,k] = ind2sub(ndim,aind);

brain_vol = single([i;j;k]');

FWHM = 4;
c = FWHM/(2*sqrt(2*log(2)));
gsp=2*c^2;

gaussianTotal = zeros(size(brain_vol,1),1,'single');
for ie = 1:size(data,1)
    elec_map_ind = pdist2(brain_vol,data(ie,1:3))<=2*FWHM;
    elec_map_vol = brain_vol(elec_map_ind,:);
    b_z = abs(elec_map_vol(:,3) - data(ie,3));
    b_y = abs(elec_map_vol(:,2) - data(ie,2));
    b_x=abs(elec_map_vol(:,1) - data(ie,1));
    gaussian = data(ie,4)*exp((-(b_x.^2+b_z.^2+b_y.^2))/gsp); %gaussian
    %     calculate the sum of electrode density step by step
    gaussianTotal(elec_map_ind) = gaussianTotal(elec_map_ind) + gaussian;
end
gaussianTotal(abs(gaussianTotal)<0.1) = 0;
data_interpolated(aind) = gaussianTotal;

end

function [data_interpolated,Ic] = ccep_interpolate_ln(ndim,data,dmax,n,method)
%n: Output image spatial resolution [mm]
% [x,y,z] = meshgrid(bb(1,1):n:bb(2,1),...
%     bb(1,2):n:bb(2,2),...
%     bb(1,3):n:bb(2,3));

data_interpolated = nan*zeros(ndim,'single');
aind = 1:prod(ndim);
[x,y,z] = ind2sub(ndim,aind);

brain_vol = single([x;y;z]');

% n1=length(1:ndim(1));
% n2=length(1:ndim(2));
% n3=length(1:ndim(3));
Ic = [];
dmax = dmax^2;
for ie = 1:size(data,1)
    d = (brain_vol(:,1)-data(ie,1)).^2+(brain_vol(:,2)-data(ie,2)).^2+(brain_vol(:,3)-data(ie,3)).^2;
    Ic = [Ic; find(d<dmax)];
end
Ic = unique(Ic);
vc = double(brain_vol(Ic,:));
% data_interpolated = NaN*zeros(n2,n1,n3); % changed by Liang 2017/04/23
data_interpolated(Ic) = griddata(data(:,1),data(:,2),data(:,3),data(:,4),vc(:,1),vc(:,2),vc(:,3),method);

end

