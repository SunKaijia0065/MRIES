function ccep_comp_batch(varargin)


if isempty(varargin)
    datapath = uigetdir('','Please choose the subject directory');
else
    datapath=varargin{1};
end



files = dir([datapath filesep 'data' filesep 'ccep*.mat']);



if exist([datapath filesep 'subjectinfo.txt'],'file')
    fid = fopen([datapath filesep 'subjectinfo.txt']);
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

prompt = {'Number of Contact per electrde [16 14 10 16 14 12 10]',...
    'Label for electrode{A;B;A^;B^}' ...
    'Bad contacts [100 131]'...
    'Number of stimulation pulse [50]'...
    'Sampling rate [2000]'};
dlg_title = 'Electrode Information';
num_lines = 1;
default_answer = { num2str(chan_array), cell2mat(chan_label),'[]','[50]','[2000]'};

if isempty(varargin)
    answer = inputdlg(prompt, dlg_title, num_lines, default_answer);
else
    answer=varargin{2};
end





numelec = str2num(answer{1});
% chan_range_info = answer{2};
% brace1 = strfind(chan_range_info, '{');
% brace2 = strfind(chan_range_info, '}');
% semicolon = strfind(chan_range_info, ';');
% cutinfo = [brace1 semicolon brace2];
% for i = 1:length(cutinfo)-1
%     chan_per_elec{i} = str2num(chan_range_info(cutinfo(i)+1:cutinfo(i+1)-1));
% end

for i = 1:length(numelec)
    chan_per_elec{i} = 1:numelec(i);
end
total_elecs = 1:sum(numelec);
label_range_info = answer{2};
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

badelec = str2num(answer{3});%from the input
numstim = str2num(answer{4});
Fs = str2num(answer{5});

% chan_name=cell(sum(chan_per_elec),1);
chan_name = cell(length(cell2mat(chan_per_elec)),1);
n=1;
for i = 1:length(chan_ID)
    for j = chan_per_elec{i}
        chan_name(n) = {strcat(chan_ID{i},num2str(j))};
        n=n+1;
    end
    %%
end
trigelec = sum(numelec)+1;
if ~exist([ datapath filesep 'stimulationdata'],'dir')
    mkdir([datapath filesep 'stimulationdata']);
end
% sum_chan = [0 cumsum(chan_per_elec)];
win = [0.4 0.6];
sort_files = natsortfiles({files.name})';
for j =1:length(sort_files)
    [~, filename, ext] = fileparts(sort_files{j});
    indx = strfind(filename,'_');
    elec1 = filename(indx(end-1)+1:indx(end)-1);
    elec2 = filename(indx(end)+1:end);
    fprintf('Stim on %s %s electrodes: ',elec1,elec2);
%     if elec1==103
%         numstim=49;
%     end

     ccep_comp_cal(datapath,[str2num(elec1) str2num(elec2)],numstim,Fs,win,badelec,trigelec,numelec,chan_name);
%     if elec1==103
%         numstim=50;
%     end
end
save([datapath filesep 'subj_elec_info.mat'],'numelec','chan_per_elec','total_elecs','win','badelec','trigelec','chan_name','Fs');
cd([datapath filesep 'stimulationdata'])