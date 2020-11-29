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
replacetime=5; %artif
data(badelec,:)=0;


%% low pass filter
lower_bound_bt = 0.1;
total_elecs = 1:size(data,1)-1;
good_elecs=setdiff(total_elecs,[badelec trigelec]);
for ielec = good_elecs
%     fprintf('filtering elec %d \n', ielec);

    [data(ielec,:)]=ft_preproc_highpassfilter(data(ielec,:),Fs,lower_bound_bt,4,[],'twopass');
end

data = data - ones(size(data,1),1)*mean(data);


%% epoching
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
    end
    
end


%% delete bad trials
ccep = [];

elec_bad_trials = cell(size(data,1),1);
total_elecs = 1:size(data,1);
val_elecs=setdiff(total_elecs,trigelec);


%% response detection

ccep_HF = cell(length(val_elecs),1);
ccep_arti=cell(length(val_elecs),1);
for ielec = val_elecs
    
    epoch_data = datamat(data(ielec,:)',locs_sec,Fs,win);
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
    %% ccep data matrix
    ccep_mat{ielec,1} = epoch_data;
    ccep_arti{ielec,1} = epoch_data_HF;
    
    ccep = [ccep mean(epoch_data,2)];
    
end

%% ccep response matrix
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
    
    envelopAll(ielec,:)=envelopMean;
    [snrP(ielec,1),latency_snr(ielec,1),latency_snrPeak(ielec,1)]=SNR(sig_test_ccep_envelopO,Fs,win);
    end
    
    %     if mean(sig_test_ccep((win(1)+0.002)*Fs:end,:),2)>0 | mean(sig_test_ccep((win(1)+0.002)*Fs:end,:),2)<0
    %         title (['REJECTED: elec' num2str(ielec)],'color','b')
    %         continue
    %     end
    ccep_elec = ccep(:,ielec);


    slidingtime = 0.010;
    %% Z score,significant point by point
    win0=0.1;
    
    sig_test_ccep=sig_test_ccep(ceil(size(sig_test_ccep,1)*(1-(win0+win(2))/(win(1)+win(2)))):end,:);
    bs_mean = mean(sig_test_ccep(1:win0*Fs-slidingtime*Fs,:),2);
    bs_std = std(bs_mean,0,1);
    m_ccep = mean(sig_test_ccep,2);
    z = (m_ccep-mean(bs_mean))./bs_std;
    pval = 1-cdf('Normal',z,0,1);%
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
                val = [val;(-1)*peaks];    
                val_indx = [val_indx;p(locs)];
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
            k(1:val_indx) = 0;
            rightindx = find(k, 1 ) - val_indx;
            if isempty(rightindx)
                rightindx = inf;
            end
            fwhm_P1(ielec,1) = leftindx + rightindx;
            
            
        end
    end
    stat.h(:,ielec) = h;    
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
    


    % permutation test by JNM 2019 method
    logN=log(shMean);
    logR=log(varAll);
    snrP=1-normcdf(logR,mean(logN),std(logN));
    
    
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
    
    