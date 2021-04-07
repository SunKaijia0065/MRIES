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
fontsize = 18;
%intergret
screensize = get(0,'MonitorPosition');
ui.mainwindow = figure('Visible', 'On', 'Name', 'MRIES', 'NumberTitle', 'Off', 'MenuBar', 'none', 'Position', ...
    [1/2*screensize(3)-600,1/2*screensize(4)-300,400,900], 'CloseRequestFcn', @mainwindow_CloseRequestFcn);



ui.subinfo = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.1 0.95 0.8 0.03], ...
    'String','Subject & Parameter','FontSize',fontsize,'Fontname', fontname);
ui.info = uipanel('Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.05 0.75 0.9 0.2]);

ui.dirpath = uicontrol('Style', 'Edit', 'Parent', ui.info, 'Units', 'Normalized', 'Position', [0.1 0.55 0.8 0.35], ...
    'Callback',@dirpath_Callback,'FontSize',fontsize,'Fontname', fontname,'max',2);
strButton = '<html><center>Subject<br> Information  ';
ui.dirbrowse = uicontrol('Style', 'Pushbutton', 'Parent', ui.info, 'Units', 'Normalized', 'Position', [0.1 0.1 0.38 0.4], ...
    'String', strButton, 'Callback', @dirbrowse_Callback,'FontSize',fontsize,'Fontname', fontname);
strButton = '<html><center>Parameter<br> Settings  ';
ui.setbrowse = uicontrol('Style', 'Pushbutton', 'Parent', ui.info, 'Units', 'Normalized', 'Position', [0.52 0.1 0.38 0.4], ...
    'String', strButton, 'Callback', @setbrowse_Callback,'FontSize',fontsize,'Fontname', fontname);

ui.interPipeling = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.1 0.7 0.8 0.03], ...
    'String','Intergreted Propressing Pipeline','FontSize',fontsize,'Fontname', fontname);
ui.interpanel = uipanel('Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.05 0.59 0.9 0.11]);



ui.pipeline = uicontrol('Style', 'Pushbutton', 'Parent', ui.interpanel, 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.8], ...
    'String', ' Data Processing Pipeline', 'Callback', @pipeline_Callback,'FontSize',fontsize,'Fontname', fontname);



%split
ui.splitFunction = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.1 0.54 0.8 0.03], ...
    'String','Separate Propressing Function','FontSize',fontsize,'Fontname', fontname);
ui.splitpanel = uipanel('Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.05 0.19 0.9 0.35]);

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
subinfo.setting  =  [];
subinfo.answer   =  [];
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

prompt = {'Number of Contact per electrde: eg. [16 14 10 16 14 12 10]',...
    'Label for electrode : eg.{A;B;A^;B^}' ...
    'Bad contacts:  eg. [100 131]'...
    'Number of stimulation pulse: eg. [50]'...
    'Sampling rate (Hz):   eg. 2000 '...
    'Stimulation interval (ms):  eg.1'};
dlg_title = 'Electrode Information';
num_lines = 1;
default_answer = { num2str(chan_array), cell2mat(chan_label),'[]','50','2000','1'};
answer = inputdlg(prompt, dlg_title, num_lines, default_answer);
subinfo.answer=answer;
fprintf(['Patient Information:','\n'])
% subinfo.answer
fprintf(['Number of contact per electrde :  ',subinfo.answer{1},'\n'])
fprintf(['Label for electrode :             ',subinfo.answer{2},'\n'])
fprintf(['Bad contacts :                    ',subinfo.answer{3},'\n'])
fprintf(['Number of stimulation pulse :     ',subinfo.answer{4},'\n'])
fprintf(['Sampling rate :                   ',subinfo.answer{5},'\n'])
fprintf(['Stimulation interval :            ',subinfo.answer{6},'\n'])


fprintf(['Completed patient initialization.\n'])
end

function setbrowse_Callback(hObject, eventdata)
    global  subinfo 
    prompt = {'Maker Channel : ''DC10'' ',...
        'High-pass filter to remove linear drift (Hz): eg. 0.1 '...
        'Baseline Time (s):  eg. 0.2 ' ,...
        'n-times std to remove the bad trial: eg. 3',...
        'Z threshold of significant N1/N2: eg. 6',...
        'Time range of significant N1/N2 (s) : eg. [0.007 0.05] s',...
        'Low response duration range (s): eg.[0.005 0.05]s',...
        'Artifict time to remove (ms) : 5',...
        'High frecquence response range (Hz): eg. : [70 170]',...
        'RMS time range (s): eg.[0.007 0.3]'};
    dlg_title = ' setting';
    num_lines = 1;
    default_answer = { 'DC10','0.1','0.2','3','6','[0.007 0.05]','[0.005 0.05]','5','[70 170]','[0.007 0.3]'};
    answer = inputdlg(prompt, dlg_title, num_lines, default_answer);
    subinfo.markerChannel = answer{1};
    subinfo.highPass      = str2num(answer{2});
    subinfo.baseTime      = str2num(answer{3});
    subinfo.zremove       = str2num(answer{4});
    subinfo.zthre         = str2num(answer{5});
    subinfo.sigN1Time     = str2num(answer{6});
    subinfo.duraRange     = str2num(answer{7});
    subinfo.replacetime   = str2num(answer{8});
    subinfo.highBand      = str2num(answer{9});  
    subinfo.rmsrange      = str2num(answer{10});

    
fprintf(['Parameter Settings:','\n'])
% subinfo.answer
fprintf(['Maker Channel :                             ',answer{1},      '\n'])
fprintf(['High-pass filter to remove linear drift :   ',answer{2},'Hz', '\n'])
fprintf(['Baseline Time :                             ',answer{3},'s', '\n'])
fprintf(['n-times std to remove the bad trial :       ',answer{4},      '\n'])
fprintf(['Z threshold of significant N1/N2:           ',answer{5},      '\n'])
fprintf(['Time range of significant N1/N2:            ',answer{6} ,'s', '\n'])
fprintf(['Low response duration range:                ',answer{7} ,'s', '\n'])
fprintf(['Artifict time to remove:                    ',answer{8}, 'ms','\n'])
fprintf(['High frecquence response range :            ',answer{9} ,'Hz' ,'\n'])
fprintf(['RMS time range :                            ',answer{10} ,'s' ,'\n'])

end

function pipeline_Callback(hObject, eventdata)

global  subinfo 
fprintf(['--------------------------Converting and Epoching------------ ----------------------\n'])
ccep_edf2mat_ALL(subinfo.mainpath)
fprintf(['---------------------------Response Detection-------------------------------------\n'])
ccep_comp_batch()
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

fprintf(['--------------------------Response Detection------------ ----------------------\n'])
ccep_comp_batch()
fprintf(['Response Detection   done!\n'])
end
function connectivity_Callback(hObject, eventdata)
global  subinfo 
fprintf(['--------------------------Connectivity Calculation------------ ----------------------\n'])
ccep_connectivity_matrix(subinfo.mainpath)
fprintf(['Connectivity Calculation   done!\n'])
end
function visualization_Callback(hObject, eventdata)
global  subinfo 
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


