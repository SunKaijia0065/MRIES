%
%
%
%
%
%
%
%
%
%

%读取各个刺激数据，计算所有的连接矩阵
% name_all={'chengsongtao','guobin','liuzhiming','lixiangju','liyihan','pangruiyang','wangjieyang','wangmingyue','wangyaobei','yangyonghao','zhangningyu','zhuwangyuan'};
name_all={'Zhangningyu','Zhuwangyuan',};
for i=1:length(name_all)
    clearvars -except name_all i
    name=name_all{i};
    name
    method_all={'RMS_rare'};
%     method_all={'Amplitude','z_Amplitude','Latency','RMS','RMS_rare','RMS_6times'};
    %% 提取个人信息
    subject.name = [];
    subject.sex = 'male';
    subject.age = [];
    subject.sphere = 'lh';
    subject.numchan = [];
    subject.num_elec = [];
    subject.chan_array = [];
    subject.Fs = 2000;
    subject.method = 'RMS';
    subject.peak = [];
    subject.connmat = [];
    subject.T = 0;
    subject.roi_radius = 3;
    
    mainpath = ['/media/sunkj/SUN1/20190422ccep_new/' name];
    
    subject.mainpath = mainpath;
    
    coordinates = load([subject.mainpath filesep 'brain3D' filesep 'autocoordinates.mat']);
    subject.coordinates = coordinates.savecoors(:,end-2:end);
    subject.dist = dist(subject.coordinates');
    
    fprintf('Loading relevant files... \n')
    
    fprintf('Loading subject''s infomation...\n')
    if exist([mainpath filesep 'subjectinfo.txt'],'file')
        fid = fopen([mainpath filesep 'subjectinfo.txt']);
        tempinfo = textscan(fid,'%s');
        fclose(fid);
        tempinfo = tempinfo{1};
        subject.name = tempinfo{2};
        
        subject.sex = tempinfo{4};
        
        subject.age = tempinfo{6};
        
        subject.sphere = tempinfo{8};
        
        subject.num_elec = str2num(tempinfo{10});
        subject.chan_array = cellfun(@str2num,tempinfo(12:12+str2num(tempinfo{10})-1));
        subject.chan_array = subject.chan_array';
        subject.cum_chan_array = cumsum(subject.chan_array);
        subject.numchan = sum(subject.chan_array);
        
        temp_name = cell2mat(tempinfo(12+str2num(tempinfo{10})+1:end));
        
        
        cutinfo = strfind(temp_name, ';');
        for i = 1:length(cutinfo)
            if i == 1
                subject.chan_name{i} = temp_name(1:cutinfo(i)-1);
            else
                
                subject.chan_name{i} = temp_name(cutinfo(i-1)+1:cutinfo(i)-1);
            end
        end
        
        subject.conn_char = cell(subject.cum_chan_array(end),1);
        
        n = 0;
        for ien = 1:subject.num_elec
            for j = 1:subject.chan_array(ien)
                n = n + 1;
                subject.conn_char{n} = [subject.chan_name{ien} num2str(j)];
            end
        end
    else
        fprintf('Subject Infomation file not found,Please load manually.\n')
    end
    
    % subject.stand_coors = importdata([mainpath filesep 'brain3D' filesep 'MNI152_coordinates_ras.txt']);
    %
    %
    fprintf('Loading stimulation data... \n')
    if exist([mainpath filesep 'stimulationdata'],'file')
        stimfiles = dir([mainpath filesep 'data' filesep 'ccep*.mat']);
        
        elec_pairs= [];
        for j = 1:length(stimfiles)
            if isempty(strfind(stimfiles(j).name,'bad'))
                [~, filename, ~] = fileparts(stimfiles(j).name);
                indx = strfind(filename,'_');
                elec1 = filename(indx(end-1)+1:indx(end)-1);
                elec2 = filename(indx(end)+1:end);
                elec_pairs = [elec_pairs;[str2num(elec1)  str2num(elec2)]];
            end
        end
        [~,idx] = sort(elec_pairs(:,1));
        elec_pairs = elec_pairs(idx,:);
        subject.elec_pairs = elec_pairs;
        
        cum_chan_per_elec = [0 subject.cum_chan_array];
        
        for j=1:length(elec_pairs)
            elec_pairs_label{j,1} = [subject.conn_char{elec_pairs(j,1)} ,'-', subject.conn_char{elec_pairs(j,2)}];
        end
        
        
        
    else
        warndlg('Stimulation data missing,please check...')
    end
    %
    %     fprintf('Loading individual pacellation file... \n')
    
    % if exist([subject.mainpath filesep 'mri' filesep 'aparc+aseg.nii'],'file')
    %     mri_info = MRIread([subject.mainpath filesep 'mri' filesep 'aparc+aseg.nii']); % input argument
    %     subject.aparc_vol = int32(mri_info.vol); % use int to save space and increase speed
    % elseif exist([subject.mainpath filesep 'mri' filesep 'aparc+aseg.mgz'],'file')
    %     mri_info = MRIread([subject.mainpath filesep 'mri' filesep 'aparc+aseg.mgz']); % input argument
    %     subject.aparc_vol = int32(mri_info.vol); % use int to save space and increase speed
    % else
    %     fprintf('Failed to load parcellation file,corresponding info may unable to display.\n')
    % end
    %
    % if exist([subject.mainpath filesep 'mri' filesep 'brain.nii'],'file')
    %     vol=MRIread([subject.mainpath filesep 'mri' filesep 'brain.nii']);
    % elseif exist([subject.mainpath filesep 'mri' filesep 'brain.mgz'],'file')
    %     vol=MRIread([subject.mainpath filesep 'mri' filesep 'brain.mgz']);
    % end
    % subject.vox2ras = vol.vox2ras;
    % subject.ras2vox = inv(subject.vox2ras);
    % subject.vox2ras_tkr = [-1,0,0,128;0,0,1,-128;0,-1,0,128;0,0,0,1];
    % subject.ras2vox_tkr = inv(subject.vox2ras_tkr);
    % fprintf('Generating standard glass brain... \r\n')
    % [ui,subject] = glass_brain_standard(ui,subject);
    % hold(ui.mainaxisc)
    % hold(ui.mainaxiss)
    % hold(ui.mainaxisa)
    %
    % ui.stim_elec_selectc = scatter(0,0,1,'filled','Parent',ui.mainaxisc,'MarkerEdgeColor','w','linewidth',1);
    % ui.stim_elec_selects = scatter(0,0,1,'filled','Parent',ui.mainaxiss,'MarkerEdgeColor','w','linewidth',1);
    % ui.stim_elec_selecta = scatter(0,0,1,'filled','Parent',ui.mainaxisa,'MarkerEdgeColor','w','linewidth',1);
    %
    % % ui.stim_elec_otherc = scatter(0,0,1,'filled','Parent',ui.mainaxisc,'MarkerEdgeColor','w','linewidth',1);
    % ui.stim_elec_others = scatter(0,0,1,'filled','Parent',ui.mainaxiss,'MarkerEdgeColor','w','linewidth',1);
    % ui.stim_elec_othera = scatter(0,0,1,'filled','Parent',ui.mainaxisa,'MarkerEdgeColor','w','linewidth',1);
    
    %     if exist([subject.mainpath filesep 'mri' filesep 'aparc+aseg.nii'],'file')
    %         mri_info = load_nii([subject.mainpath filesep 'mri' filesep 'aparc+aseg.nii']); % input argument
    %         subject.aparc_vol = mri_info.img; % use int to save space and increase speed
    %     else
    %         fprintf('Failed to load parcellation file,corresponding info may unable to display.\n')
    %     end
    
    % setup orthogonal view for T1 image
    %     if ~exist([subject.mainpath filesep 'brain3D' filesep 'T1_reslice.nii'],'file')
    %          reslice_nii([subject.mainpath filesep 'brain3D' filesep 'T1.nii'],[subject.mainpath filesep 'brain3D' filesep 'T1_reslice.nii'],[],0);
    %     end
    
    
    
    fprintf('Complete. \n')
    
    
    for j=1:length(method_all)
        subject.method=method_all{j};
        calculate_connectivity_matrix(subject);
        j
    end
    
end