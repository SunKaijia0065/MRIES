function connectivity = calculate_connectivity_matrix(subinfo)

datapath = [subinfo.mainpath filesep 'stimulationdata'];
chan_per_elec = subinfo.chan_array;
cum_chan_per_elec = cumsum(chan_per_elec);
num_chan_per_elec = sort([1 cum_chan_per_elec cum_chan_per_elec(1:end-1)+1]);
connectivity = zeros(sum(chan_per_elec));
% time range of RMS (sec) after the stimuluation
cceprange = [0.007 0.3];
% time windows for the epoched data
win = [0.1 0.5];
% sampling rate
Fs = 2000;
cceprange_point = cceprange*Fs;
win_point = win*Fs;
indx = win_point(1)+(cceprange_point(1):cceprange_point(2));
%%
% cont_per_chan = [];
% for i = 1:length(subinfo.chan_array)
%     cont_per_chan{i} = 1:subinfo.chan_array(i);
% end
% n=1;
% chan_name = cell(length(cell2mat(cont_per_chan)),1);
%
% for i = 1:length(subinfo.chan_name)
%     for j = cont_per_chan{i}
%         chan_name(n) = {strcat(subinfo.chan_name{i},num2str(j))};
%         n=n+1;
%     end
% end
chan_name = subinfo.conn_char;
%%
for i = 1:sum(chan_per_elec)-1
    i
    stimelec = [i,i+1];
    %     indxj = strfind(chan_name{stimelec(1)},';');
    elec1 = chan_name{stimelec(1)};
    %     elec1(indxj) = [];
    %     indxj = strfind(chan_name{stimelec(2)},';');
    elec2 = chan_name{stimelec(2)};
    %     elec2(indxj) = [];
    filename = [datapath filesep 'ccep_elec_' elec1 '_' elec2, '_All.mat'];
    if exist(filename,'file')
        load(filename);
        
        switch subinfo.method %没有响应的为nan ，没有数据的为0
            case 'SNR'
                
                snrPlog=-log(snrP);
                connectivity(i,:)=snrPlog;
            case 'L_SNR'
                connectivity(i,:)=latency_snr;
            case 'L_SNRPeak'
                connectivity(i,:)=latency_snrPeak;
            case 'Amplitude'
                % Amplitude
                %
                %% evoked timing range
                fwhm_dur_all(fwhm_dur_all<= 100) = 1;
                peak_final(fwhm_dur_all ~=1) = nan;
                
                connectivity(i,:) = abs(peak_final);
            case 'z_Amplitude'
                %% evoked timing range
                fwhm_dur_all(fwhm_dur_all<= 100) = 1;
                z_final(fwhm_dur_all ~=1) = nan;
                
                connectivity(i,:) = abs(z_final);
            case 'Latency'
                %% evoked timing range
                fwhm_dur_all(fwhm_dur_all<= 100) = 1;
                latency_final(fwhm_dur_all ~=1) = nan;
                connectivity(i,:) = latency_final;
            case 'RMS'
                % RMS
                h = zeros(size(all_ccep,1),size(all_ccep,2));
                if stimelec(1)    == sum(chan_per_elec)-1
                    h(1:size(stat.h,1),1:sum(chan_per_elec)-2) = stat.h(:,1:end-2);
                else
                    h(1:size(stat.h,1),:) = stat.h;
                end
                all_ccep_mask = all_ccep.*h;
                %                 plot((1:length(all_ccep_mask(:,1)))/2-100,all_ccep_mask(:,1))
                %                 hold on
                %                 plot(indx/2-100,zeros(length(indx),1),'r')
                %                 close all
                rms = sqrt(sum(all_ccep_mask(indx,:).^2)./length(indx));
                %% evoked timing range
                fwhm_dur_all(fwhm_dur_all<= 100) = 1;
                
                rms(fwhm_dur_all ~=1) = nan;
                connectivity(i,:) = rms;
            case 'RMS_rare'
                % RMS
                h = zeros(size(all_ccep,1),size(all_ccep,2));
                if stimelec(1)    == sum(chan_per_elec)-1
                    h(1:size(stat.h,1),1:sum(chan_per_elec)-2) = stat.h(:,1:end-2);
                else
                    h(1:size(stat.h,1),:) = stat.h;
                end
                all_ccep_mask = all_ccep;
                %                 plot((1:length(all_ccep_mask(:,1)))/2-100,all_ccep_mask(:,1))
                %                 hold on
                %                 plot(indx/2-100,zeros(length(indx),1),'r')
                %                 close all
                rms = sqrt(sum(all_ccep_mask(indx,:).^2)./length(indx));
                %% evoked timing range

                connectivity(i,:) = rms;
            case 'RMS_6times'
                % RMS
                h = zeros(size(all_ccep,1),size(all_ccep,2));
                if stimelec(1)    == sum(chan_per_elec)-1
                    h(1:size(stat.h,1),1:sum(chan_per_elec)-2) = stat.h(:,1:end-2);
                else
                    h(1:size(stat.h,1),:) = stat.h;
                end
                   h_indx=[zeros(indx(1)-1,stimonset_elec-1);h(indx,:);zeros(1200-indx(end),stimonset_elec-1)]; 
                
                %                 plot((1:length(all_ccep_mask(:,1)))/2-100,all_ccep_mask(:,1))
                %                 hold on
                %                 plot(indx/2-100,zeros(length(indx),1),'r')
                %                 close all
                rms = sqrt(sum((all_ccep.*h_indx).^2)./sum(h_indx));
                rms(isnan(rms))=0;
                %% evoked timing range
                fwhm_dur_all(fwhm_dur_all<= 100) = 1;
                
                rms(fwhm_dur_all ~=1) = nan;
                connectivity(i,:) = rms;
        end
        artifact_range_mm = 5;
        artifact_elecs_mm = setdiff(find(subinfo.dist(stimelec(1),:) < artifact_range_mm |  subinfo.dist(stimelec(2),:) < artifact_range_mm), stimelec);
        artifact_elecs = artifact_elecs_mm;
        connectivity(i,artifact_elecs) = 0;
    else
        connectivity(i,:) = 0;
    end
end
% supply the last electrode
connectivity(sum(chan_per_elec),:) = 0;

%the same contact
index=find(eye(size(connectivity))==1);

index=[index',index'-1];

connectivity(index(index~=0))=0;
% remove the response of the same electrode after stimuluating a pair contacts

% for j = 1:2:length(num_chan_per_elec)
%     connectivity(num_chan_per_elec(j):num_chan_per_elec(j+1),num_chan_per_elec(j):num_chan_per_elec(j+1))=0;
% end

% % % transpose to display
% % connectivity = connectivity';

if ~exist([ subinfo.mainpath filesep 'brain3D'],'dir')
    mkdir([subinfo.mainpath filesep 'brain3D']);
end
save([subinfo.mainpath filesep 'brain3D' filesep 'conn_matrix_N1+P1_' subinfo.method '.mat'],'connectivity');

