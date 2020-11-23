function ccep_connectivity_matrix(varargin)


if isempty(varargin)
    mainpath = uigetdir('','Please choose the subject directory');
else
    mainpath=varargin{1};
end

coordinates = load([mainpath filesep 'brain3D' filesep 'autocoordinates.mat']);
coordinates = coordinates.savecoors(:,end-2:end);
contactdist = dist(coordinates');
    
    if exist([mainpath filesep 'subjectinfo.txt'],'file')
        fid = fopen([mainpath filesep 'subjectinfo.txt']);
        tempinfo = textscan(fid,'%s');
        fclose(fid);
        tempinfo = tempinfo{1};
        name = tempinfo{2};
        
        sex = tempinfo{4};
        
        age = tempinfo{6};
        
        sphere = tempinfo{8};
        
        num_elec = str2num(tempinfo{10});
        chan_array = cellfun(@str2num,tempinfo(12:12+str2num(tempinfo{10})-1));
        chan_array = chan_array';
        cum_chan_array = cumsum(chan_array);
        

        
        temp_name = cell2mat(tempinfo(12+str2num(tempinfo{10})+1:end));
        
        
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
        
    else
        fprintf('Subject Infomation file not found,Please load manually.\n')        
    end
    
    
    
    
    
    
    chan_name=conn_char;
datapath = [mainpath filesep 'stimulationdata'];
chan_per_elec = chan_array;
cum_chan_per_elec = cumsum(chan_per_elec);
num_chan_per_elec = sort([1 cum_chan_per_elec cum_chan_per_elec(1:end-1)+1]);


connectivitySNR= zeros(sum(chan_per_elec));
connectivitySNR_L =zeros(sum(chan_per_elec));
connectivityA = zeros(sum(chan_per_elec));
connectivityZA = zeros(sum(chan_per_elec));
connectivityL = zeros(sum(chan_per_elec));
connectivityRMS = zeros(sum(chan_per_elec));
connectivityRMSR= zeros(sum(chan_per_elec));
        
        
% time range of RMS (sec) after the stimuluation
cceprange = [0.007 0.3];
% time windows for the epoched data
win = [0.1 0.5];
% sampling rate
Fs = 2000;
cceprange_point = cceprange*Fs;
win_point = win*Fs;
indx = win_point(1)+(cceprange_point(1):cceprange_point(2));




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
                connectivitySNR(i,:)=snrPlog;
                
            %% latency_peak
                connectivitySNR_L(i,:)=latency_snrPeak;
            %% amplitude peak
                fwhm_dur_all(fwhm_dur_all<= 100) = 1;
                peak_final(fwhm_dur_all ~=1) = nan;
                
                connectivityA(i,:) = abs(peak_final);
            %% z_amplitude
                fwhm_dur_all(fwhm_dur_all<= 100) = 1;
                z_final(fwhm_dur_all ~=1) = nan;
                
                connectivityZA(i,:) = abs(z_final);
            %% latency
                %% evoked timing range
                fwhm_dur_all(fwhm_dur_all<= 100) = 1;
                latency_final(fwhm_dur_all ~=1) = nan;
                connectivityL(i,:) = latency_final;
            %% RMS
                % RMS
%                 h = zeros(size(all_ccep,1),size(all_ccep,2));
%                 if stimelec(1)    == sum(chan_per_elec)-1
%                     h(1:size(stat.h,1),1:sum(chan_per_elec)-2) = stat.h(:,1:end-2);
%                 else
%                     h(1:size(stat.h,1),:) = stat.h;
%                 end
%                 all_ccep_mask = all_ccep.*h;

                %                 plot((1:length(all_ccep_mask(:,1)))/2-100,all_ccep_mask(:,1))
                %                 hold on
                %                 plot(indx/2-100,zeros(length(indx),1),'r')
                %                 close all
                rms = sqrt(sum(all_ccep(indx,:).^2)./length(indx));
                %% evoked timing range
                fwhm_dur_all(fwhm_dur_all<= 100) = 1;
                
                rms(fwhm_dur_all ~=1) = nan;
                connectivityRMS(i,:) = rms;
             %% RMS rare
                % RMS
                all_ccep_mask = all_ccep;
                %                 plot((1:length(all_ccep_mask(:,1)))/2-100,all_ccep_mask(:,1))
                %                 hold on
                %                 plot(indx/2-100,zeros(length(indx),1),'r')
                %                 close all
                rms = sqrt(sum(all_ccep_mask(indx,:).^2)./length(indx));
                connectivityRMSR(i,:) = rms;
               
        artifact_range_mm = 5;
        artifact_elecs_mm = setdiff(find(contactdist(stimelec(1),:) < artifact_range_mm |  contactdist(stimelec(2),:) < artifact_range_mm), stimelec);
        artifact_elecs = artifact_elecs_mm;
        
        connectivitySNR(i,artifact_elecs) = 0;
        connectivitySNR_L(i,artifact_elecs) = 0;
        connectivityA(i,artifact_elecs) = 0;
        connectivityZA(i,artifact_elecs) = 0;
        connectivityL(i,artifact_elecs) = 0;
        connectivityRMS(i,artifact_elecs) = 0;
        connectivityRMSR(i,artifact_elecs) = 0;
        
        
        
    else
        connectivitySNR(i,:) = 0;
        connectivitySNR_L(i,:) = 0;
        connectivityA(i,:) = 0;
        connectivityZA(i,:) = 0;
        connectivityL(i,:) = 0;
        connectivityRMS(i,:) = 0;
        connectivityRMSR(i,:) = 0;
    end 
    
    
    

end

%the same contact

connectivtiys={connectivitySNR,connectivitySNR_L,connectivityA,connectivityZA,...
    connectivityL,connectivityRMS,connectivityRMSR};
connectivityNames={'SNR','L_SNRPeak','Amplitude','z_Amplitude','Latency','RMS','RMS_rare'};

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


