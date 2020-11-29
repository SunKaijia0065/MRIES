% Convert data from edf format to mat
% Writen by Yunxian Bai at 20171220
% changed by Wang Liang at 20180129
% updata by Kaijia Sun
% data structure:
% # hospital:
%       # subject:
%             # edf folder (storing edf files, which can be splited by two electrodes)
%
function ccep_edf2mat_bh(varargin)
Fs = 2000;
if isempty(varargin)
    [filename,pathname] = uigetfile('*.*',  'Please select the edf file to be converted');
else
    filename = varargin{1};
    pathname = varargin{2};
end
fprintf('%s\n',['The selected file is ' filename]);
cd(pathname)
[hdr, data] = edfRead(filename);
Annotation = data(end,:);
label = hdr.label;
label = cellfun(@(x) strrep(x,'EEG ', ''), label, 'UniformOutput', false);
label = cellfun(@(x) strrep(x,'POL ', ''), label, 'UniformOutput', false);
label = cellfun(@(x) strrep(x,'-Ref', ''), label, 'UniformOutput', false);
delete_channel = {'E', 'EDF Annotations','POL'};
data(ismember(label,delete_channel),:)= [];
label(ismember(label,delete_channel))= [];

DC_channel = find(ismember(label,{'DC10'}));
data = data([1:DC_channel-1,DC_channel+1:end, DC_channel],:);
label = label([1:DC_channel-1,DC_channel+1:end, DC_channel]);

elec_not_ABC = find(ismember(label,{'ECG','EMG1','EMG2','EMG3','EMG4','EMG5'}));
data(elec_not_ABC,:) = [];
label(elec_not_ABC) = [];

if exist(['..' filesep 'subjectinfo.txt'],'file')
    fid = fopen(['..' filesep 'subjectinfo.txt']);
    for n = 1:4
        tline = fgetl(fid);
        
    end
    tempinfo = textscan(fid,'%s');
    fclose(fid);
    
    num_elec = str2num(tempinfo{1}{2});
    chan_array = cellfun(@str2num,tempinfo{1}(4:4+num_elec-1))';
    chan_label = tempinfo{1}(4+num_elec+1:end)';
else
    fprintf('Subject Infomation file not found. Pleae create it first manually.\n')
end

label_range_info = chan_label{:};
brace1 = strfind(label_range_info, '{');
brace2 = strfind(label_range_info, '}');
if any(strfind(label_range_info, ';'))
    semicolon = strfind(label_range_info, ';');
elseif any(strfind(label_range_info, ' '))
    semicolon = strfind(label_range_info, ' ');
end
cutinfo = [brace1 semicolon brace2];
for i = 1:length(cutinfo)
    if i == 1
        chan_ID{i} = label_range_info(1:cutinfo(i)-1);
    else
        
        chan_ID{i} = label_range_info(cutinfo(i-1)+1:cutinfo(i)-1);
    end
end
% index of A' B'...A B...G
indx = [];
for i = 1:length(chan_ID)
    order  = [];
    k= [] ;
    for j = 1:length(label)
        
        Num=str2double(regexp(label{j},'\d*','match')');
        n = strfind(label{j},num2str(Num));
        if strcmp(label{j}(1:n-1),chan_ID{i})
            k = [k;str2num(label{j}(length(chan_ID{i})+1:end))];
            order = [order;j];
        end
        
    end
    [~,m] = sort(k);
    indx = [indx;order(m)];
end
index = [indx;size(data,1)];


default_chan_per_elec = chan_array;

label = label(index);
whole_data = data(index,:);

clear data;


prompt = {'Number of contact per electrde: eg. [16 14 10 16 14 14 14 8 10 10]',...
    'Label of electrode, eg. [1 2]',...
    'Contact of electrode, eg. {[1:10, 13:15]; [2:8, 10:13]}' };
dlg_title = 'parameters';
num_lines = 1;
default_answer = { ['[' num2str(default_chan_per_elec) ']'], '[]', '{[]}'};

if isempty(varargin)
    answer = inputdlg(prompt, dlg_title, num_lines, default_answer);
else
    answer = { ['[' num2str(default_chan_per_elec) ']'], ['',varargin{3},''], ['{',varargin{4},'}']}';
    
end

chan_per_elec = str2num(answer{1});
elec_range = str2num(answer{2});

chan_range_info = answer{3};
brace1 = strfind(chan_range_info, '{');
brace2 = strfind(chan_range_info, '}');
semicolon = strfind(chan_range_info, ';');
cutinfo = [brace1 semicolon brace2];
for i = 1:length(cutinfo)-1
    chan_range{i} = str2num(chan_range_info(cutinfo(i)+1:cutinfo(i+1)-1));
end


cum_chan_per_elec = [0 cumsum(chan_per_elec)];

stim_elec = [];
for i = 1:length(elec_range)
    s = [chan_range{i}' chan_range{i}'+1];
    %     if chan_range{i}(end) == chan_per_elec(elec_range(i))
    %         s(length(chan_range{i}),:)=[];
    %     end
    stim_elec = [stim_elec;cum_chan_per_elec(elec_range(i))+s];
end

if ~exist('../data')
    mkdir('../data');
end
cd('../data')

mx = max(Annotation);
mn = min(Annotation);
if abs(mx) > abs(mn)
    [~,locs] = findpeaks(double(Annotation),'MINPEAKDISTANCE',0.5*Fs,'MinPeakHeight',abs(mx)*0.8);
else
    [~,locs] = findpeaks((-1).*double(Annotation),'MINPEAKDISTANCE',0.5*Fs,'MinPeakHeight',abs(mn)*0.8);
end
difference = diff(locs);
ind = find(difference>Fs*50&difference<Fs*54); % You many need to change 50 or 54 if the data length is not in the range (sec): 50-54

if size(stim_elec,1) ~= length(ind)
    error('Please Check the edf_match.txt and .edf data\n')
end
if length(ind)==1
    fprintf('%s\n',['Generating file ccep_elec_' num2str(stim_elec(1,1)) '_' num2str(stim_elec(1,2)) '.mat'])
    data = whole_data;
    save(['ccep_elec_' num2str(stim_elec(1,1)) '_' num2str(stim_elec(1,2))], 'data');
else
    for i = 1:length(ind)
        fprintf('%s\n',['Generating file ccep_elec_' num2str(stim_elec(i,1)) '_' num2str(stim_elec(i,2)) '.mat'])
        if locs(ind(i))<1000
            data = whole_data(:,1:locs(ind(i)+1)+4000);
        elseif size(whole_data,2)<locs(ind(i)+1)+4000
            data = whole_data(:,locs(ind(i))-1000:end);
        else
            data = whole_data(:,locs(ind(i))-1000:locs(ind(i)+1)+4000);
        end
        % output the valid pair stimuluation
        save(['ccep_elec_' num2str(stim_elec(i,1)) '_' num2str(stim_elec(i,2))], 'data');
    end
end
%% if you already install EEGLAB, you can check one of the data using the following command
% eegplot(data(1:end,:),'srate',Fs,'winlength',10,'dispchans',20,'spacing',3200,'position',[30 30 1500 1200]);