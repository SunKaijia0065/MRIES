%Computing the connectivity matrix from response detection file
% updated by Kaijia Sun
function ccep_connectivity_matrix(varargin)
%patient information
global subinfo
if isempty(varargin)
    mainpath = uigetdir('','Please choose the subject directory');
else
    mainpath=varargin{1};
end


coordinates = load([mainpath filesep 'brain3D' filesep 'autocoordinates.mat']);
coordinates = coordinates.savecoors(:,end-2:end);
contactdist = dist(coordinates');


if ~isempty(subinfo)
    chan_array = str2num(subinfo.answer{1});
    num_elec = length(chan_array);
    temp_name = subinfo.answer{2};
    
elseif exist([mainpath filesep 'subjectinfo.txt'],'file')
    fid = fopen([mainpath filesep 'subjectinfo.txt']);
    tempinfo = textscan(fid,'%s');
    fclose(fid);
    tempinfo = tempinfo{1};
    num_elec = str2num(tempinfo{10});
    chan_array = (cellfun(@str2num,tempinfo(12:12+str2num(tempinfo{10})-1)))';
    temp_name = cell2mat(tempinfo(12+str2num(tempinfo{10})+1:end));
    
else
    fprintf('Subject Infomation file not found,Please load manually.\n')
end

    cum_chan_array = cumsum(chan_array);    

    
    
cutinfo = strfind(temp_name, ';');
for i = 1:length(cutinfo)
    if i == 1
        chan_name{i} = temp_name(1:cutinfo(i)-1);
    else

        chan_name{i} = temp_name(cutinfo(i-1)+1:cutinfo(i)-1);
    end
end

conn_char = cell(cum_chan_array(end),1);

n = 0;
for ien = 1:num_elec
    for j = 1:chan_array(ien)
        n = n + 1;
        conn_char{n} = [chan_name{ien} num2str(j)];
    end
end
    



chan_name=conn_char;
datapath = [mainpath filesep 'stimulationdata'];
chan_per_elec = chan_array;
cum_chan_per_elec = cumsum(chan_per_elec);
num_chan_per_elec = sort([1 cum_chan_per_elec cum_chan_per_elec(1:end-1)+1]);


connectivity_HF_str = zeros(sum(chan_per_elec));
connectivity_HF_lat = zeros(sum(chan_per_elec));
connectivity_LF_amp = zeros(sum(chan_per_elec));
connectivity_LF_Zamp= zeros(sum(chan_per_elec));
connectivity_LF_lat = zeros(sum(chan_per_elec));
connectivity_LF_sRMS= zeros(sum(chan_per_elec));
connectivity_LF_RMS = zeros(sum(chan_per_elec));


% time range of RMS (sec) after the stimuluation
rmsrange = subinfo.rmsrange;%RMS time range
% time windows for the epoched data
baseTime=subinfo.baseTime;%baseTime
% sampling rate
Fs = str2num(subinfo.answer{5});%read
fwhm_durRange=subinfo.duraRange*Fs;%

rmsrange_point = rmsrange*Fs;
basePoint = baseTime*Fs;
indx = basePoint+(rmsrange_point(1):rmsrange_point(2));



for i = 1:sum(chan_per_elec)
    
    stimelec = [i,i+1];
    %     indxj = strfind(chan_name{stimelec(1)},';');
    elec1 = chan_name{stimelec(1)};
    %     elec1(indxj) = [];
    %     indxj = strfind(chan_name{stimelec(2)},';');
    if stimelec(2)>sum(chan_per_elec)
        elec2=elec1;
    else
        elec2 = chan_name{stimelec(2)};
    end
    fprintf(['loading stimulation   ',elec1,'-',elec2,'\n']);
    %     elec2(indxj) = [];
    filename = [datapath filesep 'ccep_elec_' elec1 '_' elec2, '_All.mat'];
    if exist(filename,'file')
        load(filename);
        
        
        %% snr
        
        snrPlog=-log(snrP);
        connectivity_HF_str(i,:)=snrPlog;
        
        %% latency_peak
        connectivity_HF_lat(i,:)=latency_snrPeak;
        %% amplitude peak
        fwhm_dur_all((fwhm_dur_all<= fwhm_durRange(2))&(fwhm_dur_all>= fwhm_durRange(1))) = 1;
        peak_final(fwhm_dur_all ~=1) = nan;
        
        connectivity_LF_amp(i,:) = abs(peak_final);
        %% z_amplitude
        fwhm_dur_all((fwhm_dur_all<= fwhm_durRange(2))&(fwhm_dur_all>= fwhm_durRange(1))) = 1;
        z_final(fwhm_dur_all ~=1) = nan;
        
        connectivity_LF_Zamp(i,:) = abs(z_final);
        %% latency
        %% evoked timing range
        fwhm_dur_all((fwhm_dur_all<= fwhm_durRange(2))&(fwhm_dur_all>= fwhm_durRange(1))) = 1;
        latency_final(fwhm_dur_all ~=1) = nan;
        connectivity_LF_lat(i,:) = latency_final;
        %% RMS
        rms = sqrt(sum(all_ccep(indx,:).^2)./length(indx));
        %% evoked timing range
        fwhm_dur_all((fwhm_dur_all<= fwhm_durRange(2))&(fwhm_dur_all>= fwhm_durRange(1))) = 1;
        
        rms(fwhm_dur_all ~=1) = nan;
        connectivity_LF_sRMS(i,:) = rms;
        %% RMS rare
        % RMS
        all_ccep_mask = all_ccep;
        rms = sqrt(sum(all_ccep_mask(indx,:).^2)./length(indx));
        connectivity_LF_RMS(i,:) = rms;
        
        %% distance to delete
        artifact_range_mm = 5;
        artifact_elecs_mm = setdiff(find(contactdist(stimelec(1),:) < artifact_range_mm |  contactdist(stimelec(2),:) < artifact_range_mm), stimelec);
        artifact_elecs = artifact_elecs_mm;
        
        connectivity_HF_str(i,artifact_elecs) = 0;
        connectivity_HF_lat(i,artifact_elecs) = 0;
        connectivity_LF_amp(i,artifact_elecs) = 0;
        connectivity_LF_Zamp(i,artifact_elecs) = 0;
        connectivity_LF_lat(i,artifact_elecs) = 0;
        connectivity_LF_sRMS(i,artifact_elecs) = 0;
        connectivity_LF_RMS(i,artifact_elecs) = 0;
        
        
        
    else
        connectivity_HF_str(i,:) = 0;
        connectivity_HF_lat(i,:) = 0;
        connectivity_LF_amp(i,:) = 0;
        connectivity_LF_Zamp(i,:) = 0;
        connectivity_LF_lat(i,:) = 0;
        connectivity_LF_sRMS(i,:) = 0;
        connectivity_LF_RMS(i,:) = 0;
    end
    
    
    
    
end

%the same contact

connectivtiys={connectivity_HF_str,connectivity_HF_lat,connectivity_LF_amp,connectivity_LF_Zamp,...
    connectivity_LF_lat,connectivity_LF_sRMS,connectivity_LF_RMS};
connectivityNames={'HF_str','HF_lat','LF_amp','LF_Zamp','LF_lat','LF_sRMS','LF_RMS'};

for i = 1:length(connectivtiys)
    connectivity=connectivtiys{i};
    method=connectivityNames{i};
    connectivity=sti_con_delete(connectivity);
    
    if ~exist([ mainpath filesep 'Matrix'],'dir')
        mkdir([mainpath filesep 'Matrix']);
    end
    save([mainpath filesep 'Matrix' filesep 'conn_matrix_' method '.mat'],'connectivity');
    
end


end


function connectivity=sti_con_delete(connectivity)
index=find(eye(size(connectivity))==1);

index=[index',index'-1];

connectivity(index(index~=0))=0;
end


