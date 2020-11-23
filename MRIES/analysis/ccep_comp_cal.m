function [z_N1,peak_N1, latency_N1] = ccep_comp_cal(datapath,stimelec,numstim,Fs,win,badelec,trigelec,chan_per_elec,chan_name)
% Compute cortico-cortical evoked potential by averaging over multiple
% waves around the stimulated electrodes
%
%Inputs:
% datapath: full path for storing the segmented data
% stimelec: paired elctrodes exposed to the brain stimulus (eg. [35 36])
% numstim: number of brain stimuli given (eg. 49)
% Fs: sampling rate (eg. 2000)
% win: time range (sec) (eg. [0.1 0.5], which means left 0.1s (not negative value) and right 0.5s around stim onset)
% badelec: bad channels (eg. [127 128])
% trigelec: electrode with trigger signal
% chan_per_elec: number of contacts for each electrode
% chan_name: channel name
%
%Requirment:
% 1) should install EEGlab before running this code
% 2) data file was segmented to ONLY include one paired electrodes stimululated per time (eg., ccep_elec_35_36.mat with "data" as a variable)
%
% Author: Liang Wang
% Date: 2015-05-10
%
% Changed on 2018-01-01 Liang Wang



cd([datapath filesep 'data'])
filename = ['ccep_elec_' num2str(stimelec(1)) '_' num2str(stimelec(2)), '.mat'];
if nargin < 5
    error('Please add more input parameters \n');
end
if nargin < 7
    trigelec = [];
end

badelec = [badelec stimelec];


if exist(filename,'file')
    load(filename)
    
else
    error('This file does not exist \n')
end
fprintf('%s \n', filename);

data = double(data); % row: electrodes; column: time courses
% decompression to 1/5
% for i = 1:size(data,1)
%     datanew(i,:)=resample(data(i,:),1,5);
% end
% data = datanew;
% clear datanew;
% Fs = Fs/5;
replacetime=5;
data(badelec,:)=0;

%% band-stop filter to remove power line noises
% data = data;
% lower_bound_bt = [ 48 97 146 196];
% upper_bound_bt = [ 52 103 154 200];
%
% for ielec = 1:size(data,1)-1
% %     fprintf('filtering elec %d \n', ielec);
% %     [data(ielec,:)]=band_stop(data(ielec,:), Fs, lower_bound_bt, upper_bound_bt, 1);
%     [data(ielec,:)]=band_stop_new(data(ielec,:), Fs, size(data,2), 50, 50, 3, 20);
% end

%% low pass filter
lower_bound_bt = 0.1;
% upper_bound_bt = [ 200];
total_elecs = 1:size(data,1)-1;
good_elecs=setdiff(total_elecs,[badelec trigelec]);
for ielec = good_elecs
%     fprintf('filtering elec %d \n', ielec);

    [data(ielec,:)]=ft_preproc_highpassfilter(data(ielec,:),Fs,lower_bound_bt,4,[],'twopass');
end

data = data - ones(size(data,1),1)*mean(data);


%% band pass filter
%     lower_bound_bt = [ 0.1 ];
%     upper_bound_bt = [ 200 ];
% 
%     for ielec = 1:size(data,1)
% 
%         [data(ielec,:)]=napl_band_pass(data(ielec,:), Fs, lower_bound_bt, upper_bound_bt, 1);
%     end

%% re-reference data
% total_elecs = 1:size(data,1);
% all_elecs=setdiff(total_elecs,[badelec trigelec]);% defines good electrodes for the file-wide CAR
% CAR_bank=sum(data(all_elecs,:),1)/length(all_elecs); %calculates CAR
%
% gdat_clean=data;
% for e = all_elecs
%     gdat_clean(e,:) = data(e,:) - CAR_bank;%subtracts bank CAR from each electrode in the bank
% end
%
% data = gdat_clean;
% clear gdat_clean;

% check raw data
% eegplot(data(1:end,:),'srate',Fs,'winlength',10,'dispchans',20,'spacing',3200,'title',filename(1:end-4),'position',[30 30 1500 1200]);


peakflg = 1;
while peakflg
    if isempty( trigelec)
        stimonset_elec = input('Select an electrode with clear stim effect []:');
    else
        stimonset_elec = trigelec;
        if size(data) > trigelec
            data(trigelec+1:end,:)=[];
        end
    end
    
    if ismember(stimonset_elec,stimelec)
        peakflg = 1;
        warning('Please re-select another electrode with more clear stim effect');
        trigelec= trigelec+1;
        if trigelec>sum(chan_per_elec)
            trigelec=trigelec-sum(chan_per_elec);
        end
        fprintf('Trigger eletrode is %s \n',num2str(trigelec));
        continue
    end
    mx = max(data(stimonset_elec,:));
    mn = min(data(stimonset_elec,:));
    if abs(mx) > abs(mn)
        [pks,locs] = findpeaks(data(stimonset_elec,:),'MINPEAKDISTANCE',0.9*Fs,'SORTSTR','descend','NPEAKS',numstim);
        if length(locs) == numstim
            peakflg = 0;
        else
            peakflg = 1;
            warning('Please re-select another electrode with more clear stim effect');
            trigelec= trigelec+1;
            if trigelec>sum(chan_per_elec)
                trigelec=trigelec-sum(chan_per_elec);
            end
            fprintf('Trigger eletrode is %s \n',num2str(trigelec));
            continue
        end
    else
        [pks,locs] = findpeaks((-1).*data(stimonset_elec,:),'MINPEAKDISTANCE',0.9*Fs,'SORTSTR','descend','NPEAKS',numstim);
        if length(locs) == numstim
            peakflg = 0;
        else
            peakflg = 1;
            warning('Please re-select another electrode with more clear stim effect');
            trigelec= trigelec+1;
            if trigelec>sum(chan_per_elec)
                trigelec=trigelec-sum(chan_per_elec);
            end
            fprintf('Trigger eletrode is %s \n',num2str(trigelec));
            continue
        end
    end
    
    
    %%
    [locs, indx]=sort(locs);
    
    pks = pks(indx);
    for i =1:length(locs) 
        beforeIndex=locs(i)-10:locs(i);
        seg=data(stimonset_elec,beforeIndex);
        segIndex=find(seg<min(seg)+0.001*(max(seg)-min(seg)));
        locs(i)=locs(i)-10-1+segIndex(end);
    end
    %show the start time
%     plot(data(stimonset_elec,:))
%     hold on
%     scatter(locs,zeros(length(locs),1)+0)
    
    
    
    interval = diff(locs);
    locs_test = interval > 0.9*Fs & interval < 1.25*Fs;
    if any(~locs_test) || locs(end) > (size(data,2)-win(2)*Fs)
        peakflg = 1;
        warning('Please re-check the sorting process')
        trigelec= trigelec+1;
        if trigelec>sum(chan_per_elec)
            trigelec=trigelec-sum(chan_per_elec);
        end
        fprintf('Trigger eletrode is %s \n',num2str(trigelec));
        
    else
        peakflg = 0;
%         fprintf('%d peaks \n', length(locs));
    end
    
end

%% low pass filter
% lower_bound_bt = [ 1 ];
% upper_bound_bt = [ 40 ];
% total_elecs = 1:size(data,1)-1;
% good_elecs=setdiff(total_elecs,[badelec trigelec]);
% locs_data = zeros(size(data));
% for ielec = good_elecs
% %     fprintf('filtering elec %d \n', ielec);
%     for i = 1:length(locs)
%         locs_data(ielec,locs(i)-4:locs(i)+4) = data(ielec,locs(i)-4:locs(i)+4);
%         step = (data(ielec,locs(i)+4)-data(ielec,locs(i)-4))/8;
%         if step == 0
%             temp = ones(1,9)*data(ielec,locs(i)-4);
%         else
%             temp = data(ielec,locs(i)-4):step:data(ielec,locs(i)+4);
%         end
%         data(ielec,(locs(i)-4):(locs(i)+4)) = temp;
%     end
%     [data(ielec,:)]=band_pass(data(ielec,:), Fs, lower_bound_bt, upper_bound_bt, 1);
% end
% data = data+locs_data;

% eegplot(data,'srate',Fs,'winlength',10,'dispchans',20,'spacing',300,'title',filename(1:end-4),'position',[30 30 1500 1200])
locs_sec = locs./Fs;
badstim = [];
locs_sec(badstim)=[];


%% bipolar ���㷽ʽ
% bipolardata = ecog_bipolarize(data,chan_per_elec);
% data_new = [bipolardata;data(end,:)];
% data = data_new;
%%
ccep = [];

elec_bad_trials = cell(size(data,1),1);
num_stim = length(locs_sec);
n1 = round(num_stim/2);
n1_indx = 1:n1;%& randperm(num_stim,n1);
n2_indx = setdiff(1:num_stim,n1_indx);

stimonsetEpcoh=datamat(data(stimonset_elec,:)',locs_sec(n1_indx),Fs,win);
% plot(stimonsetEpcoh)

total_elecs = 1:size(data,1);
val_elecs=setdiff(total_elecs,trigelec);

good_elecs=setdiff(total_elecs,[badelec trigelec]);
%% split the trials in tow groups to test the replication
%ccep_1
ccep1_mat = cell(length(val_elecs),1);
ccep1_arti = cell(length(val_elecs),1);
ccep_HF = cell(length(val_elecs),1);
ccep_arti=cell(length(val_elecs),1);
for ielec = val_elecs
    
    epoch_data = datamat(data(ielec,:)',locs_sec,Fs,win);
    
        %vis the signal
%     time=((1:size(epoch_data(:,1),1))/Fs-0.1)*1000;
%     ind=find((time<20)&(time>-10));
%     plot(time(ind),epoch_data(ind,1))

    
    baseline = mean(epoch_data((win(1)-0.2)*Fs:win(1)*Fs-1,:));
    epoch_data = epoch_data - ones(size(epoch_data,1),1)*baseline;
    %% the high-frecquence
    epoch_data_HF=replace_artifict(epoch_data,Fs,replacetime,win);
%     for i = 1:size(epoch_data_HF,2)
%         epoch_data_HFadd(:,i)=bandpass(epoch_data_HF(:,i),[70 170 ],Fs);
%     end
    %% remove bad trials
    badtrials=remwave(epoch_data);
    %         figure,plot(epoch_data);
    %         figure,plot(epoch_data(:,largewave));
    elec_bad_trials{ielec} = badtrials;
    epoch_data(:,badtrials) = [];
    epoch_data_HF(:,badtrials) = [];
    %%
    ccep_mat{ielec,1} = epoch_data;
    ccep_arti{ielec,1} = epoch_data_HF;
    
    ccep = [ccep mean(epoch_data,2)];
    
end
% ccep1(round((win(1)-0.002)*Fs):round((win(1)+0.002)*Fs),:) = 0; % remove the stim effect (4 ms)


%% plot ccep1


% nrow = length(chan_per_elec);
% ncol = max(chan_per_elec);
% col_ind = 1:ncol:nrow*ncol;
%
% e_trans = [];
% for ielec = 1:nrow
%     e_trans = [e_trans col_ind(ielec)+chan_per_elec(ielec)-1:-1:col_ind(ielec)];
% end
% ScSz = get( groot, 'Screensize' );
%
% pos_grid =  [ScSz(1)+30 ScSz(2)+30  ScSz(3)-30  ScSz(4)-30];
% grid1 = figure('Name',['ccep-elec-' num2str(stimelec(1)) '-' num2str(stimelec(2))]);
% set(grid1,'Position', ScSz,'color',[1 1 1],'MenuBar','figure');
% gridPos = DivideScreen(nrow,ncol,ScSz,40,60,pos_grid);
peak_N1 = zeros(sum(chan_per_elec),1); latency_N1 = zeros(sum(chan_per_elec),1); z_N1= zeros(sum(chan_per_elec),1);
peak_P1 = zeros(sum(chan_per_elec),1); latency_P1 = zeros(sum(chan_per_elec),1); z_P1= zeros(sum(chan_per_elec),1);
snrP=zeros(sum(chan_per_elec),1);
latency_snr=zeros(sum(chan_per_elec),1);
latency_snrPeak=zeros(sum(chan_per_elec),1);
fwhm_P1 = zeros(sum(chan_per_elec),1); fwhm_N1 = zeros(sum(chan_per_elec),1);

stat.p = [];
stat.h = [];
stat.z = [];
min_amplitude = 30; % uV
for ielec = val_elecs
    %     axes('Position',gridPos{e_trans(ielec)});
%      fprintf(['elec #' num2str(ielec) '\n']  );
    
    
    %     if ismember(ielec,stimelec)
    %         title(['STIM: ' chan_name{ielec} ] ,'color','r')
    %         continue
    %     end
    %     if ismember(ielec,badelec)
    %         title (['BAD:' chan_name{ielec}],'color','b')
    %         continue
    %     end
    sig_test_ccep = ccep_mat{ielec};
    if isempty(sig_test_ccep)
        ccep_HF{ielec}=nan;
        ccep_arti{ielec}=nan;
        snrP(ielec,1)=nan;
        latency_snr(ielec,1)=nan;
        latency_snrPeak(ielec,1)=nan;
    else
    
    sig_test_ccep_arti=ccep_arti{ielec};
    % remove long-lasting artifacts by removing the average signal
    sig_test_ccep_arti=sig_test_ccep_arti-repmat(mean(sig_test_ccep_arti,2),1,size(sig_test_ccep_arti,2));
    ccep_arti{ielec}=sig_test_ccep_arti;
    sigOne=reshape(sig_test_ccep_arti,1,numel(sig_test_ccep_arti));
    %     sigBand=bandpass(sigOne,[70 170 ],Fs);%signal
    sigBand=ft_preproc_bandpassfilter(sigOne,Fs,[70 170],8,[],'twopass');
    test=reshape(sigBand,size(sig_test_ccep_arti));
    
    
    enveloponly=abs(hilbert(sigBand));%abs of signal
    %     enveloponly=ft_preproc_lowpassfilter(enveloponly,Fs,20,[],[],'twopass');
    %     sig_test_ccep_HF=reshape(sigBand,size(sig_test_ccep_arti));
    sig_test_ccep_HF=reshape(enveloponly,size(sig_test_ccep_arti));
    ccep_HF{ielec}=sig_test_ccep_HF;
    sig_test_ccep_envelopO=reshape(enveloponly,size(sig_test_ccep_arti));
    bs = mean(sig_test_ccep_envelopO((win(1)-0.2)*Fs:win(1)*Fs-1,:));
    sig_test_ccep_envelopO = sig_test_ccep_envelopO - ones(size(sig_test_ccep_envelopO,1),1)*bs;
    
    envelopMean=mean(sig_test_ccep_envelopO,2);
    % plot the envelop
%     plot(sig_test_ccep_HF,'g')
%     hold on
%     plot(envelopMean,'r')
    
    envelopAll(ielec,:)=envelopMean;
    [snrP(ielec,1),latency_snr(ielec,1),latency_snrPeak(ielec,1)]=SNR(sig_test_ccep_envelopO,Fs,win);
    end
    
    %     if mean(sig_test_ccep((win(1)+0.002)*Fs:end,:),2)>0 | mean(sig_test_ccep((win(1)+0.002)*Fs:end,:),2)<0
    %         title (['REJECTED: elec' num2str(ielec)],'color','b')
    %         continue
    %     end
    ccep_elec = ccep(:,ielec);

    %     if abs(min(min(ccep_elec))) > max(max(ccep_elec))
    %         ylim([min(min(ccep_elec))-10 -min(min(ccep_elec))+10]);
    %     else
    %         ylim([-max(max(ccep_elec))-10 max(max(ccep_elec))+10]);
    %     end
    
    %     if abs(min(min(ccep_elec))) < 300 && abs(max(max(ccep_elec))) < 300
    %         ylim([-300 300]);
    %     end
    %     ylim([-1000 1000]);
    %     ylimit = get(gca,'Ylim');
    %     plot([0 0],[ylimit(1) ylimit(2)],'-k');
    
    slidingtime = 0.010;
    % test the significant area
    %% sign test
    
    
    
    %     bs = mean(sig_test_ccep(1:win(1)*Fs,:));
    %     p = zeros(size(sig_test_ccep,1)-slidingtime*Fs/2,1);
    %     h = zeros(size(sig_test_ccep,1)-slidingtime*Fs/2,1);
    %     z = zeros(size(sig_test_ccep,1)-slidingtime*Fs,1);
    %     compN = length(round((win(1)+0.002)*Fs) : size(sig_test_ccep,1)-slidingtime*Fs/2);
    %     for j = round((win(1)+0.002)*Fs) : size(sig_test_ccep,1)-slidingtime*Fs/2
    %         slideing_ccep = mean(sig_test_ccep(j-slidingtime*Fs/2+1 : j + slidingtime*Fs/2,:));
    %         [p(j,1),h(j,1), s] = signrank(slideing_ccep,bs,'alpha', 0.05/compN);
    %         if ~isfield(s,'zval')
    %             s.zval = 0;
    %         end
    %         z(j,1) = s.zval;
    %     end
    %     h(1:round((win(1)+0.002)*Fs),:) = 0;
    
    %% Z score
    win0=0.1;
    
    sig_test_ccep=sig_test_ccep(ceil(size(sig_test_ccep,1)*(1-(win0+win(2))/(win(1)+win(2)))):end,:);
    bs_mean = mean(sig_test_ccep(1:win0*Fs-slidingtime*Fs,:),2);%�������?
    bs_std = std(bs_mean,0,1);%��׼��
    m_ccep = mean(sig_test_ccep,2);
    z = (m_ccep-mean(bs_mean))./bs_std;
    pval = 1-cdf('Normal',z,0,1);%����һ���ֲ�������p_value����
    h = abs(z) > 6; % 6 SD
    h(1:round((win0+0.002)*Fs),:) = 0;
    %% Find connected components in the significant indices (at least lasting 10ms)
    CC = bwconncomp(h,4);
    for icluster = 1:CC.NumObjects
        if length(CC.PixelIdxList{icluster}) <= slidingtime*Fs
            h(CC.PixelIdxList{icluster}) = 0;
        end
    end
    stat.p(:,ielec) = pval;
    stat.h(:,ielec) = h;
    stat.z(:,ielec) = z;
    %% find the earlier N1 latency and amplitude
    CC_new = bwconncomp(h,4);
    %     m_ccep = mean(sig_test_ccep,2);
    all_ccep(:,ielec) = m_ccep;
    max_N1_indx = win0*Fs + 0.05*Fs; %% N1 is less than 0.05 sec
    min_N1_indx = win0*Fs + 0.007*Fs; %% N1 is larger than 0.007 sec
    %     ncount = 1;
    val = [];val_indx = [];
    for icluster = 1:CC_new.NumObjects
        p = CC_new.PixelIdxList{icluster};
        q = p(p<max_N1_indx & p>min_N1_indx);
        %         [peaks locs]=min(m_ccep(p));
        [peak_i loc_i] = findpeaks((-1)*m_ccep(p));
        [peaks indx_i] = max(peak_i);
        locs = loc_i(indx_i);
        if mean(m_ccep(p)) < 0 && ~isempty(q)
            if max(peaks) > min_amplitude % a peak does exceed 20 uV
                %             [mn k]=min(m_ccep(q));
                %             val(ncount,1) = mn;
                %             val_indx(ncount,1)= q(k);
                val = [val;(-1)*peaks];    %��ȡֵ��λ��
                val_indx = [val_indx;p(locs)];
                %             ncount = ncount + 1;
            else
                h(p) = 0;
            end
        end
    end
    if ~isempty(val) && ~isempty(val_indx)
        val = val(val_indx<=max_N1_indx & val_indx>=min_N1_indx);
        val_indx = val_indx(val_indx<=max_N1_indx & val_indx>=min_N1_indx);
        
        if ~isempty(val) && ~isempty(val_indx)
            if length(val)>1
                [Y,I] = max(abs(val));
                val = val(I);
                val_indx = val_indx(I);
                
            end
            peak_N1(ielec,1) = val(1);
            latency_N1(ielec,1) = val_indx(1)/Fs-win0;
            z_N1(ielec,1)=z(val_indx(1));
            %% find the full width at half-maximum
            k=m_ccep > val(1)/2;
            m = k;
            m(val_indx:end) = 0;
            leftindx = abs(find(m, 1, 'last' ) - val_indx);
            k(1:val_indx) = 0;
            rightindx = find(k, 1 ) - val_indx;
            
            if isempty(rightindx)
                rightindx = inf;
            end
            fwhm_N1(ielec,1) = leftindx + rightindx;
            
            %             plot(latency_N1(ielec,1),peak_N1(ielec,1),'rs','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',5)
        end
    end
    %% find the earlier P1 latency and amplitude
    max_P1_indx = win0*Fs + 0.05*Fs;
    min_P1_indx = win0*Fs + 0.007*Fs;
    val = [];val_indx = [];
    for icluster = 1:CC_new.NumObjects
        p = CC_new.PixelIdxList{icluster};
        q = p(p<max_P1_indx & p>min_P1_indx);
        %         [peaks locs] = findpeaks(m_ccep(p));
        [peak_i loc_i] = findpeaks(m_ccep(p));
        [peaks indx_i] = max(peak_i);
        locs = loc_i(indx_i);
        if mean(m_ccep(p)) > 0 && ~isempty(q)
            if max(peaks) > min_amplitude % a peak does exceed 20 uV
                val = [val;peaks];
                val_indx = [val_indx;p(locs)];
            else
                h(p) = 0;
            end
            
        end
    end
    if ~isempty(val) && ~isempty(val_indx)
%         %find good result
%         plot(z)
%         hold on
%         plot([1 1200],[6 6]);plot([1 1200],[-6 -6])
%         axis([-100 1000 -100 100])
        
        val = val(val_indx<=max_P1_indx & val_indx>=min_P1_indx);
        val_indx = val_indx(val_indx<=max_P1_indx & val_indx>=min_P1_indx);
        if ~isempty(val) && ~isempty(val_indx)
            if length(val)>1
                [Y,I] = max(abs(val));
                val = val(I);
                val_indx = val_indx(I);
                
                
            end
            peak_P1(ielec,1) = val(1);
            latency_P1(ielec,1) = val_indx(1)/Fs-win0;
            z_P1(ielec,1)=z(val_indx(1));
            %% find the full width at half-maximum
            k=m_ccep < val(1)/2;
            m = k;
            m(val_indx:end) = 0;
            leftindx = abs(find(m, 1, 'last' ) - val_indx);
            %             leftindx = min(abs(find(m_ccep(find(k)<val_indx)) - val_indx))-1;
            k(1:val_indx) = 0;
            rightindx = find(k, 1 ) - val_indx;
            if isempty(rightindx)
                rightindx = inf;
            end
            fwhm_P1(ielec,1) = leftindx + rightindx;
            
            
            %             plot(latency_P1(ielec,1),peak_P1(ielec,1),'rs','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',5)
        end
    end
    stat.h(:,ielec) = h;
    %     %%
    % %     plot(-win(1)+1/Fs:1/Fs:win(2)-slidingtime/2,h*ylimit(2)/2.*sign(z),'m');
    %     plot(-win(1)+1/Fs:1/Fs:win(2),h*ylimit(2)/2.*sign(z),'m');
    % %     area([40 800], [ylimit(2) ylimit(2)], ylimit(1), 'FaceColor', [0.5 .5 .5]);
    %
    %
    %     set(gca,'YDir','reverse')
    %     xlim([-0.1 0.5]);
    %     if h(round((win(1)+0.005)*Fs):end,:)>0 | h(round((win(1)+0.005)*Fs):end,:)<0
    %         title (['REJECTED: ' chan_name{ielec}],'color','b')
    %     else
    %         title (['' chan_name{ielec}])
    %     end
    %
    %     box off
    %     axcopy
    %
    
end
latency_N1(latency_N1==0) = nan;
latency_P1(latency_P1==0) = nan;

latency = [latency_N1 latency_P1];
peak_temp = [peak_N1 peak_P1];
z_temp=[z_N1 z_P1];

[latency_final,index] = min(latency,[],2);
for idx = 1:length(index)
    peak_final(idx,1) = peak_temp(idx,index(idx));
    z_final(idx,1) = z_temp(idx,index(idx));
end

envelopAll=envelopAll';

fwhm_dur=[fwhm_N1 fwhm_P1];
fwhm_dur(fwhm_dur==0)=nan;
[fwhm_dur_all,~] = min(fwhm_dur,[],2);
stimelec = [stimelec(1) stimelec(2)];
% rerun_rmglm = 1;

cd('../stimulationdata')
% hgsave(grid1,['ccep_elec_' chan_name{stimelec(1)} '_'  chan_name{stimelec(2)} '.fig']);
save(['ccep_elec_' chan_name{stimelec(1)} '_' chan_name{stimelec(2)} '_All.mat'],...
    'peak_N1','latency_N1','peak_P1','latency_P1','peak_final','latency_final',...
    'all_ccep','stat','stimonset_elec','stimelec','fwhm_dur','fwhm_dur_all',...
    'z_N1','z_P1','z_final','snrP','envelopAll','latency_snr','latency_snrPeak');
% save('tianjie','stimelec','stimonset_elec','-append')

fprintf('Saved...\n');
cd(datapath)

% figure('Name',['ccep1-elec-' num2str(stimelec(1)) '-' num2str(stimelec(2))])
%
% plot(-win(1)+1/Fs:1/Fs:win(2),ccep1);
% xlim([-0.1 0.5]);
% ylim([min(min(ccep1))-10 max(max(ccep1))+10]);
% if abs(min(min(ccep1))) < 300 && abs(max(max(ccep1))) < 300
%     ylim([-300 300]);
% end
% panelname = {'All good elecs showed'};
% title(panelname)
% ylabel('Amplitude (uV)'); xlabel('Time aligned to stim onset (sec)');


% ep = ccep1((win(1)+0.01)*Fs:end,:);
% ep_thr = 50;
% bad_ccep_elec=unique([find(sum(ep > ep_thr)./size(ep,1)>0.7)  find(sum(ep < -ep_thr)./size(ep,1)>0.7)]); % at least 70% exceeds the threshold values
%
% %% don't trait good elecs as bad ones
% % if the evoked potentials has less value than the threshold, which means
% % there exists kind of oscilliations for this electrode (good one)
% q = [];
% for p = 1:length(bad_ccep_elec)
%     if any(ep(:,bad_ccep_elec(p)) > 0)
%         if sum(ep(:,bad_ccep_elec(p)) < ep_thr)/size(ep,1)>0.2 % at least 20% exceeds the threshold values
%             q = [q p];
%         end
%     elseif any(ep(:,bad_ccep_elec(p)) < 0)
%         if sum(ep(:,bad_ccep_elec(p)) > -ep_thr)/size(ep,1)>0.2
%             q = [q p];
%         end
%     end
% end
% bad_ccep_elec(q)=[];
% ccep1(:,bad_ccep_elec) = 0;
%
% %%
% subplot(2,1,2)
% plot(-win(1)+1/Fs:1/Fs:win(2),ccep1);
% ep_elec=remwave(ccep1);
% if ~isempty(num2str(bad_ccep_elec))
%     %     figname = {['ccep1-elec-' num2str(stimelec(1)) '-' num2str(stimelec(2))],['elecs with ccep1:' num2str(ep_elec)],...
%     %         ['bad CCEP of elecs:' num2str(bad_ccep_elec)]};
%     panelname = {['bad CCEP of elecs removed:' num2str(bad_ccep_elec)]};
%
% else
%     %     figname = {['ccep1-elec-' num2str(stimelec(1)) '-' num2str(stimelec(2))],['elecs with ccep1:' num2str(ep_elec)]};
%     panelname = {'All elecs are good'};
%
% end
%
%
% ylim([min(min(ccep1))-10 max(max(ccep1))+10]);
% if abs(min(min(ccep1))) < 300 && abs(max(max(ccep1))) < 300
%     ylim([-300 300]);
% end
% xlim([-0.1 0.5]);
% title(panelname)
% ylabel('Amplitude (uV)'); xlabel('Time aligned to stim onset (sec)');


function data=datamat(data,E,Fs,win)
% Inputs:
% data   (input time series as a column vector) - required
% E      (events to use as triggers) - required
% Fs     (sampling frequency of data) - required
% win    (window around triggers to use data matrix -[winl winr]) - required
%          e.g [1 1] uses a window starting 1 * Fs samples before E and
%              ending 1*Fs samples after E.
% Note that E, Fs, and win must have consistent units

NE=length(E);
nwinl=round(win(1)*Fs);
nwinr=round(win(2)*Fs);
nE=floor(E*Fs)+1;
datatmp=[];
for n=1:NE
    indx=nE(n)-nwinl:nE(n)+nwinr;
    if min(indx)<=0
        fprintf('some baseline absence!\n')
        indx(indx<=0)=1;
    end
    datatmp=[datatmp data(indx)];
end
data=datatmp;
function epoch_data=replace_artifict(epoch_data,Fs,replacetime,win)
    time=-win(1):1/Fs:win(2);
    time=time(2:end);
    indArti=find((time>=0-1e-5)&(time<=replacetime/1000+1e-5));
    indBefore=find((time>=-replacetime/1000-1e-5)&(time<=0+1e-5));
    indBefore=indBefore(end:-1:1);
    mult=1:-1/(size(indBefore,2)-1):0;
    mult=repmat(mult,size(epoch_data,2),1);
    sigBefore=epoch_data(indBefore,:).*mult';
    
    indAfter=find((time>=replacetime/1000-1e-5)&(time<=(2*replacetime)/1000+1e-5));
    indAfter=indAfter(end:-1:1);
    mult=0:1/(size(indBefore,2)-1):1;
    mult=repmat(mult,size(epoch_data,2),1);
    sigAfter=epoch_data(indAfter,:).*mult';   
    epoch_data(indArti,:)=sigBefore+sigAfter;
    


function [largewave, points] = remwave(data,stdtime)
% Inputs:
% data   (input time series as a column vector) - required
%
if nargin < 2
    stdtime = 5; % 3 standard deviation
end
low=nanmean(data,2)-stdtime*nanstd(data,0,2);
high=nanmean(data,2)+stdtime*nanstd(data,0,2);
% low=nanmedian(data,2)-stdtime*nanstd(data,0,2);
% high=nanmedian(data,2)+stdtime*nanstd(data,0,2);
outlier1 = (data < repmat(low,1,size(data,2)));
outlier2 = (data > repmat(high,1,size(data,2)));
outlier = outlier1 | outlier2;
largewave = find(sum(outlier));
points = find(sum(outlier,2)>0);
% plot_fft




function [snrP,time0,timePeak]=SNR(signal,Fs,win)
    time=-win(1):1/Fs:win(2);
    time=time(2:end);
    timeBin=[10,25,40,55,70,85,100]/1000;
    
    
    varShurff=zeros((length(timeBin)-1),1000);
%     for i = 1:(length(timeBin)-1)
%         sig=signal((time>=timeBin(i)&(time<=timeBin(i+1))));
%         snrSp(i)=var(sig);
%     end
    sig=signal((time>=timeBin(1)&(time<=timeBin(7))),:);
    sigRaw=mean(sig,2);
    varAll=var(sigRaw);
    if varAll==0
        snrP=1;
        time0=nan;
        timePeak=nan;
        return
    end
    sig=repmat(sig,1,1000);
 
    shuInd=ceil(rand(size(sig,2),1)*size(sig,1));
    sh=circshift_columns(sig,shuInd);
    shPart=mat2cell(sh,size(sig,1),zeros(1000,1)+size(signal,2));
    shAve=cell2mat(cellfun(@(x) mean(x,2),shPart,'UniformOutput',false));
    shMean=var(shAve,1);
%     for i = 1:(length(timeBin)-1)
%         num=floor(size(shAve,1)/(length(timeBin)-1));
%         varShurff(i,:)=var(shAve(i*num-num+1:i*num,:),1);
%         varSp(i)=var(sigRaw(i*num-num+1:i*num));
%     end
    


    % permutation test by JNM 2019 method
%     varShurffMean=mean(varShurff,1);
%     snrShurff=shMean./varShurffMean;
%     varMean=mean(varSp);
%     snrReal=varAll./varMean;
    
%     logN=log(snrShurff);
%     
%     logR=log(snrReal);
%     snrP=1-normcdf(logR,mean(logN),std(logN));
    
    logN=log(shMean);
    logR=log(varAll);
    snrP=1-normcdf(logR,mean(logN),std(logN));
    
    %snrP=1-sum(shMean>varAll)/length(shMean);
    
    envelop=mean(signal,2);
    baseline=envelop((time<0)&(time>-200/1000));
    baselinevar=std(baseline);
    ther=mean(baseline)+baselinevar*norminv(0.999);
    time0=find((envelop(1:end-1)>ther)&(time>0)'&(time<300/1000)');
    
    testsig=envelop((time>10/1000)&(time<100/1000));
    [peak,peakindex]=findpeaks(testsig);
    maxpeakindex=peakindex(peak==max(peak));
    timePeak=peak(peak==max(peak));
%     [timePeak,halfindex]=max(envelop((time>10/1000)&(time<100/1000)));
    
    if isempty(time0)
        time0=nan;
    else 
        time0=time(time0(1));
    end
    
    if timePeak<ther
        timePeak=nan;
    else 
        timePeak=time(maxpeakindex+find(time==0/1000)+10*Fs/1000);
    end
    
    if isempty(timePeak)
        timePeak=nan;
    end
    
    
    
    
    
    
    

% function [snrP,time0,timePeak]=SNR_new(signal,Fs,win)
%     time=-win(1):1/Fs:win(2);
%     time=time(2:end);
%     timeBin=[10,25,40,55,70,85,100]/1000;
%     sig=signal((time>=timeBin(1)&(time<=timeBin(7))),:);
%     sigRaw=mean(sig,2);
%     varAll=var(sigRaw);
%     if varAll==0
%         snrP=1;
%         time0=nan;
%         timePeak=nan;
%         return
%     end
%      for i = 1:(length(timeBin)-1)
%         num=floor(size(sigRaw,1)/(length(timeBin)-1));
%         varSp(i)=var(sigRaw(i*num-num+1:i*num));
%     end
%     varMean=mean(varSp);
%     snrReal=varAll./varMean;
%     %% randmization test
%     signal_rev = flipud(signal);
%     snrShurff=zeros(1,1000);
%  
%     signal_rev = signal_rev(0.1*Fs:0.8*Fs,:);
%     for j = 1:1000
%         shuInd=ceil(rand(size(signal_rev,2),1)*size(signal_rev,1));
%         sh=circshift_columns(signal_rev,shuInd);
%         
%         time_sh=-0.3:1/Fs:0.4;
%         time_sh=time_sh(2:end);
%         timeBin=[10,25,40,55,70,85,100]/1000;
%         sig_rand=sh((time_sh>=timeBin(1)&(time_sh<=timeBin(7))),:);
%         sigRaw_rand=mean(sig_rand,2);
%         varAll_rand=var(sigRaw_rand);
%         for i = 1:(length(timeBin)-1)
%             num=floor(size(sigRaw_rand,1)/(length(timeBin)-1));
%             varSp_rand(i)=var(sigRaw_rand(i*num-num+1:i*num));
%         end
%         varMean_rand=mean(varSp_rand);
%         snrShurff(j)=varAll_rand./varMean_rand;
%     end
%     logN=log(snrShurff);
%     
%     logR=log(snrReal);
%     snrP=1-normcdf(logR,mean(logN),std(logN));
% 
%     
%     %snrP=sum(snrShurff>snrReal)/length(snrShurff);
%     
%     envelop=mean(signal,2);
%     baseline=envelop((time<0)&(time>-200/1000));
%     baselinevar=std(baseline);
%     ther=mean(baseline)+baselinevar*norminv(0.999);
%     time0=find((envelop(1:end-1)>ther)&(time>0)'&(time<300/1000)');
%     
%     [timePeak,halfindex]=max(envelop((time>0)&(time<100/1000)));
%     
%     if isempty(time0)
%         time0=nan;
%     else 
%         time0=time(time0(1));
%     end
%     
%     if timePeak<ther
%         timePeak=nan;
%     else 
%         timePeak=time(halfindex+find(time==0));
%     end
