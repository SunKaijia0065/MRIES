%MRIES main GUI
%updated by Kaijia Sun
%
%
function MRIES

warning('off')
close all
%% Setting UI Objects and subject infomation
global ui subinfo

fontname= 'Times New Roman';
fontsize = 15;
%intergret
screensize = get(0,'MonitorPosition');
ui.mainwindow = figure('Visible', 'On', 'Name', 'MRIES', 'NumberTitle', 'Off', 'MenuBar', 'none', 'Position', ...
    [1/2*screensize(3)-600,1/2*screensize(4)-400,400,900], 'CloseRequestFcn', @mainwindow_CloseRequestFcn);



ui.interPipeling = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.1 0.95 0.8 0.03], ...
    'String','Intergreted Propressing Pipeline','FontSize',fontsize,'Fontname', fontname);
ui.interpanel = uipanel('Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.1 0.6 0.8 0.35]);

ui.dirpath = uicontrol('Style', 'Edit', 'Parent', ui.interpanel, 'Units', 'Normalized', 'Position', [0.1 0.7 0.8 0.22], ...
    'Callback',@dirpath_Callback,'FontSize',fontsize,'Fontname', fontname,'max',2);
ui.dirbrowse = uicontrol('Style', 'Pushbutton', 'Parent', ui.interpanel, 'Units', 'Normalized', 'Position', [0.1 0.4 0.8 0.22], ...
    'String', ' Subject Floder', 'Callback', @dirbrowse_Callback,'FontSize',fontsize,'Fontname', fontname);

ui.pipeline = uicontrol('Style', 'Pushbutton', 'Parent', ui.interpanel, 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.22], ...
    'String', ' Data Processing Pipeline', 'Callback', @pipeline_Callback,'FontSize',fontsize,'Fontname', fontname);



%split
ui.splitFunction = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.1 0.55 0.8 0.03], ...
    'String','Separate Propressing Function','FontSize',fontsize,'Fontname', fontname);
ui.splitpanel = uipanel('Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.1 0.2 0.8 0.35]);

ui.edf2mat = uicontrol('Style', 'Pushbutton', 'Parent', ui.splitpanel, 'Units', 'Normalized', 'Position', [0.1 0.7 0.8 0.22], ...
    'String', ' Converting and Epoching', 'Callback', @edf2mat_Callback,'FontSize',fontsize,'Fontname', fontname);

ui.responseDect = uicontrol('Style', 'Pushbutton', 'Parent', ui.splitpanel, 'Units', 'Normalized', 'Position', [0.1 0.4 0.8 0.22], ...
    'String', ' Response Detection', 'Callback', @responseDect_Callback,'FontSize',fontsize,'Fontname', fontname);

ui.connectivity = uicontrol('Style', 'Pushbutton', 'Parent', ui.splitpanel, 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.22], ...
    'String', ' Connectivity Calculation', 'Callback', @connectivity_Callback,'FontSize',fontsize,'Fontname', fontname);

% visualization

ui.view3D = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.1 0.15 0.8 0.03], ...
    'String','Visualization','FontSize',fontsize,'Fontname', fontname);
ui.view3Dpanel = uipanel('Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.1 0.05 0.8 0.1]);
ui.visualization = uicontrol('Style', 'Pushbutton', 'Parent', ui.view3Dpanel, 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.8], ...
    'String', ' Visualization Interface', 'Callback', @visualization_Callback,'FontSize',fontsize,'Fontname', fontname);



% Set the mainpath on the gui
P = mfilename('fullpath');
[pathstr,~,~] = fileparts(P);
warning off MATLAB:iofun:UnsupportedEncoding
a=get(ui.dirpath,'string');
subinfo.mainpath = pathstr;
if isempty(a)
    set(ui.dirpath,'String',pathstr)
end








end



function dirpath_Callback(hObject, eventdata)

global subinfo
mainpath = get(hObject,'String');
subinfo.mainpath = mainpath;

end

function dirbrowse_Callback(hObject, eventdata)

global ui subinfo 
fprintf(['Input the patient...\n'])
subinfo.mainpath = uigetdir('','Please choose the subject directory');
fprintf(['Patient Path is :',subinfo.mainpath,'\n'])
set(ui.dirpath,'String',[subinfo.mainpath filesep])

if exist([subinfo.mainpath filesep 'subjectinfo.txt'],'file')
    fid = fopen([subinfo.mainpath filesep 'subjectinfo.txt']);
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
answer = inputdlg(prompt, dlg_title, num_lines, default_answer);
subinfo.answer=answer;
fprintf(['Patient Information:','\n'])
% subinfo.answer
fprintf(['Number of contact per electrde :',subinfo.answer{1},'\n'])
fprintf(['Label for electrode :',subinfo.answer{2},'\n'])
fprintf(['Bad contacts :',subinfo.answer{3},'\n'])
fprintf(['Number of stimulation pulse :',subinfo.answer{4},'\n'])
fprintf(['Sampling rate :',subinfo.answer{5},'\n'])

fprintf(['Completed patient initialization.\n'])
end

function pipeline_Callback(hObject, eventdata)

global  subinfo 
fprintf(['--------------------------Converting and Epoching------------ ----------------------\n'])
ccep_edf2mat_ALL(subinfo.mainpath)
fprintf(['---------------------------Response Detection-------------------------------------\n'])
ccep_comp_batch(subinfo.mainpath,subinfo.answer)
fprintf(['---------------------------Connectivity Calculation--------------------------------------\n'])
ccep_connectivity_matrix(subinfo.mainpath)
fprintf(['---------------------------All Processing Done--------------------------------------\n'])

end

function edf2mat_Callback(hObject, eventdata)
global  subinfo 
fprintf(['--------------------------Converting and Epoching------------ ----------------------\n'])
ccep_edf2mat_ALL(subinfo.mainpath)
fprintf(['Converting and Epoching   done!\n'])

end
function responseDect_Callback(hObject, eventdata)
global  subinfo 
fprintf(['--------------------------Response Detection------------ ----------------------\n'])
ccep_comp_batch(subinfo.mainpath,subinfo.answer)
fprintf(['Response Detection   done!\n'])
end
function connectivity_Callback(hObject, eventdata)
global  subinfo 
fprintf(['--------------------------Connectivity Calculation------------ ----------------------\n'])
ccep_connectivity_matrix(subinfo.mainpath)
fprintf(['Connectivity Calculation   done!\n'])
end
function visualization_Callback(hObject, eventdata)
global  subinfo CCEP3D
MRIESviewer(subinfo.mainpath)
end

function mainwindow_CloseRequestFcn(hObject,eventdata)

global ui subinfo  
delete(hObject);

ui = []; 
subinfo = [];  
clear ui subinfo
clc
fprintf('Thanks for using MRIES . \n')
end


