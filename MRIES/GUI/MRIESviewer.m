%MRIES viewer GUI
%updated by Kaijia Sun
%

function MRIESviewer(varargin)
%
%
%
%
warning('off')
close all
%% Setting UI Objects and subject infomation
global ui subinfo

subinfo.name = [];
subinfo.sex = 'male';
subinfo.age = [];
subinfo.sphere = 'lh';
subinfo.numchan = [];
subinfo.num_elec = [];
subinfo.chan_array = [];
subinfo.Fs = 2000;
subinfo.method = 'RMS';
subinfo.peak = [];
subinfo.connmat = [];
subinfo.T = 0;
subinfo.roi_radius = 3;
subinfo.connthresh= [];
%% Main window
screensize = get(0,'MonitorPosition');
ui.mainwindow = figure('Visible', 'On', 'Name', 'MRIES Viewer', 'NumberTitle', 'Off', 'MenuBar', 'none', 'Position', ...
    [1/2*screensize(3)-600,1/2*screensize(4)-400,1000,900], 'CloseRequestFcn', @mainwindow_CloseRequestFcn, ...
    'WindowButtonDownFcn', @mainwindow_WindowButtonDownFcn, 'WindowKeyPressFcn', @mainwindow_WindowKeyPressFcn);

ui.mainaxisc = axes('Parent', ui.mainwindow, 'Units', 'Normalized', 'OuterPosition', [0 0.5 0.5 0.5], ...
    'XMinorGrid', 'On', 'YMinorGrid', 'On', 'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0]);
axiscLtext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.04 0.75 0.01 0.02], ...
    'String', 'L', 'FontSize', 12);
axiscRtext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.47 0.75 0.01 0.02], ...
    'String', 'R', 'FontSize', 12);
axiscStext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.26 0.97 0.01 0.02], ...
    'String', 'S', 'FontSize', 12);
axiscItext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.26 0.53 0.01 0.02], ...
    'String', 'I', 'FontSize', 12);

ui.mainaxiss = axes('Parent', ui.mainwindow, 'Units', 'Normalized', 'OuterPosition', [0.5 0.5 0.5 0.5], ...
    'XMinorGrid', 'On', 'YMinorGrid', 'On', 'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0]);
axissLtext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.54 0.75 0.01 0.02], ...
    'String', 'P', 'FontSize', 12);
axissRtext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.97 0.75 0.01 0.02], ...
    'String', 'A', 'FontSize', 12);
axissStext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.76 0.97 0.01 0.02], ...
    'String', 'S', 'FontSize', 12);
axissItext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.76 0.53 0.01 0.02], ...
    'String', 'I', 'FontSize', 12);

ui.mainaxisa = axes('Parent', ui.mainwindow, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 0.5], ...
    'XMinorGrid', 'On', 'YMinorGrid', 'On', 'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0]);
axisaLtext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.04 0.25 0.01 0.02], ...
    'String', 'L', 'FontSize', 12);
axisaRtext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.47 0.25 0.01 0.02], ...
    'String', 'R', 'FontSize', 12);
axisaStext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.26 0.47 0.01 0.02], ...
    'String', 'A', 'FontSize', 12);
axisaItext = uicontrol('Style', 'Text', 'Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.26 0.03 0.01 0.02], ...
    'String', 'P', 'FontSize', 12);

ui.rascoord_text = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.50 0.48 0.04 0.02], ...
    'String','RAS:','FontSize',12);
ui.rascoord = uicontrol('Style', 'Edit','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.55 0.475 0.19 0.03], ...
    'String','','FontSize',12);
ui.volcoord_text = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.50 0.51 0.04 0.02], ...
    'String','VOL:','FontSize',12);
ui.volcoord = uicontrol('Style', 'Edit','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.55 0.505 0.19 0.03], ...
    'String','','FontSize',12);
ui.underlay_text = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.75 0.475 0.07 0.03], ...
    'String','Underlay:','FontSize',12);
ui.underlay_value = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.82 0.475 0.06 0.03], ...
    'String','','FontSize',12);
ui.overlay_text = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.88 0.475 0.06 0.03], ...
    'String','Overlay:','FontSize',12);
ui.overlay_value = uicontrol('Style', 'Text','Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.94 0.475 0.04 0.03], ...
    'String','','FontSize',12);

ui.mainpanel = uipanel('Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.53 0.03 0.44 0.44]);
% ui.mainpanel = uipanel('Parent', ui.mainwindow, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.9]);

ui.dir_text = uicontrol('Style', 'Text', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.02 0.88 0.98 0.1], ...
    'String', 'Please choose the subject directory.', 'FontSize', 12);

ui.dirpath = uicontrol('Style', 'Edit', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.02 0.80 0.76 0.10], ...
    'Callback',@dirpath_Callback);

ui.dirbrowse = uicontrol('Style', 'Pushbutton', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.8 0.8 0.18 0.10], ...
    'String', ' Browse', 'Callback', @dirbrowse_Callback);

% subject infomation control UI
ui.subinfo_text = uicontrol('Style', 'Text', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.01 0.70 0.4 0.06], ...
    'String', 'Subject Infomation:', 'FontSize', 12);% 'BackgroundColor', [0.84 0.84 0.84]

ui.name_text = uicontrol('Style', 'Text', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.02 0.63 0.16 0.06], ...
    'String', 'Name', 'FontSize', 12);
ui.name = uicontrol('Style', 'Edit', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.20 0.63 0.22 0.07], ...
    'FontSize', 12, 'Callback', @name_Callback);

ui.sex_text = uicontrol('Style', 'Text', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.02 0.55 0.16 0.06], ...
    'String', 'Sex', 'FontSize', 12);
ui.sex = uicontrol('Style', 'Popupmenu', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.20 0.55 0.22 0.06], ...
    'String', {'','male','female'},'Value', 1, 'FontSize', 12, 'Callback', @sex_Callback);

ui.age_text = uicontrol('Style', 'Text', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.02 0.44 0.16 0.06], ...
    'String', 'Age', 'FontSize', 12);
ui.age = uicontrol('Style', 'Edit', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.20 0.44 0.22 0.06], ...
    'FontSize', 12, 'Callback', @age_Callback);

ui.sphere_text = uicontrol('Style', 'Text', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.02 0.36 0.16 0.06], ...
    'String', 'Sphere', 'FontSize', 12);
ui.sphere = uicontrol('Style', 'Popupmenu', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.20 0.36 0.22 0.06], ...
    'String', {'','left','right','bilateral'}, 'FontSize', 12, 'Callback', @sphere_Callback);

ui.numchan_text = uicontrol('Style', 'Text', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.02 0.25 0.16 0.06], ...
    'String', 'numelec', 'FontSize', 12);
ui.numchan = uicontrol('Style', 'Edit', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.20 0.25 0.22 0.07], ...
    'FontSize', 12, 'Callback', @numchan_Callback);

ui.chan_array_text = uicontrol('Style', 'Text', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.02 0.16 0.42 0.06], ...
    'String', 'Electrode Information:', 'FontSize', 12);
ui.chan_array = uicontrol('Style', 'Edit', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.05 0.1 0.9 0.07], ...
    'FontSize', 12, 'Callback', @chan_array_Callback);
ui.chan_ABC = uicontrol('Style', 'Edit', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.05 0.02 0.9 0.07], ...
    'FontSize', 12, 'Callback', @chan_ABC_Callback);

ui.stimelec_text = uicontrol('Style', 'Text', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.5 0.7 0.4 0.06], ...
    'String', 'Stimulate Electrode', 'FontSize', 12);

ui.stimelec = uicontrol('Style', 'Popupmenu', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.5 0.64 0.4 0.06], ...
    'String', ' ', 'FontSize', 12, 'Callback', @stimelec_Callback);



ui.circle_button = uicontrol('Style', 'PushButton', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.5 0.25 0.22 0.1], ...
    'String', 'Circle Map','Enable', 'On', 'FontSize', 12, 'Callback', @circle_button_Callback);

ui.coordinate_button = uicontrol('Style', 'PushButton', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.75 0.25 0.22 0.1], ...
    'String', 'Coordinates','Enable', 'On', 'FontSize', 12, 'Callback', @coordinate_Callback);

% Displaying Button
% ui.orthobutton = uicontrol('Style', 'Pushbutton', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.7 0.4 0.25 0.1], ...
%     'String', 'Ortho View', 'FontSize', 12, 'Callback', @orthobutton_Callback);

ui.connbutton = uicontrol('Style', 'Pushbutton', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.5 0.36 0.22 0.1], ...
    'String', 'Conn Mat', 'FontSize', 12, 'Callback', @connbutton_Callback);

ui.surfacebutton = uicontrol('Style', 'Pushbutton', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.75 0.36 0.22 0.1], ...
    'String', 'Surf View', 'Enable', 'On', 'FontSize', 12, 'Callback', @surfacebutton_Callback);

ui.method_text = uicontrol('Style', 'Text', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.49 0.55 0.4 0.06], ...
    'String', 'Response Indicator', 'FontSize', 12);

ui.method = uicontrol('Style', 'Popupmenu', 'Parent', ui.mainpanel, 'Units', 'Normalized', 'Position', [0.5 0.5 0.4 0.06], ...
    'String',{'','RMS','RMS_rare','Amplitude','Latency','SNR','L_SNR','L_SNRPeak';}, 'Value', 1, 'FontSize', 12, 'Callback', @method_Callback);

% Set the mainpath on the gui
P = mfilename('fullpath');
[pathstr,~,~] = fileparts(P);
warning off MATLAB:iofun:UnsupportedEncoding
a=get(ui.dirpath,'string');
subinfo.mainpath = pathstr;
if isempty(a)
    set(ui.dirpath,'String',pathstr)
end

% Load parcellation and segmentation info
fid = fopen([pathstr filesep 'CCEP_colorLUT.txt']);
aparc_raw = textscan(fid,'%f %s %f %f %f %f');
fclose(fid);
ui.aparc_indx = cell2mat(aparc_raw(1));
ui.aparc_aseg = char(aparc_raw{2});
fid = fopen([pathstr filesep 'cortex_abbr.txt']);
ccep_raw = textscan(fid,'%f %s %f %f %f %f');
fclose(fid);
ui.ccep_cortex = char(ccep_raw{2});

viewpos = get(ui.mainwindow,'Position');



if ~isempty(varargin)
    dirbrowse_Callback(varargin{1});
end

%% circle map window
ui.circleview = figure('Visible', 'Off', 'Name', 'Circle Map', 'NumberTitle', 'Off', 'MenuBar', 'none', 'Position', ...
    [viewpos(1)+viewpos(3),viewpos(2)+400,600,600], 'CloseRequestFcn', @view_CloseRequestFcn);
ui.circleaxes1 = axes('Parent', ui.circleview, 'Visible', 'Off', 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.8], ...
    'XLim', [-15 15], 'YLim', [-15 15], 'Nextplot', 'Add');
ui.circleaxes2 = axes('Parent', ui.circleview, 'Visible', 'Off', 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.8], ...
    'XLim', [-15 15], 'YLim', [-15 15], 'Nextplot', 'Add');

%% connectivity matrix window
ui.connview = figure('Visible', 'Off', 'Name', 'Connectivity Matrix', 'NumberTitle', 'Off', 'MenuBar', 'none', 'Position', ...
    [viewpos(1)+0.5*viewpos(3),viewpos(2),600,600], 'CloseRequestFcn', @view_CloseRequestFcn);
ui.connaxes = axes('Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.08 0.15 0.82 0.78]);
ui.conn_text_s = uicontrol('Style', 'Text', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.05 0.01 0.05 0.06], ...
    'String', 'S:', 'FontSize', 16);
ui.conn_index_s = uicontrol('Style', 'Edit', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.12 0.03 0.15 0.04], ...
    'FontSize', 12, 'Enable', 'Inactive', 'Callback', @conn_index_x_Callback);
ui.conn_text_r = uicontrol('Style', 'Text', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.35 0.01 0.05 0.06], ...
    'String', 'R:', 'FontSize', 16);
ui.conn_index_r = uicontrol('Style', 'Edit', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.42 0.03 0.15 0.04], ...
    'FontSize', 12, 'Enable', 'Inactive', 'Callback', @conn_index_y_Callback);
ui.conn_text_value = uicontrol('Style', 'Text', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.6 0.03 0.15 0.04], ...
    'String', 'Value:', 'FontSize', 16);
ui.conn_value = uicontrol('Style', 'Text', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.75 0.03 0.18 0.04], ...
    'String', 'xxx', 'FontSize', 16);
ui.conn_thresh_text = uicontrol('Style', 'Text', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.91 0.50 0.09 0.05], ...
    'String', 'Thresh', 'FontSize', 12);
ui.conn_thresh_value = uicontrol('Style', 'Text', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.92 0.48 0.06 0.03], ...
    'String', ' ', 'FontSize', 10);
ui.conn_thresh = uicontrol('Style', 'Slider', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.91 0.16 0.08 0.32], ...
    'Callback', @conn_thresh_Callback);


ui.conn_dist_text = uicontrol('Style', 'Text', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.91 0.92 0.09 0.05], ...
    'String', 'Dist', 'FontSize', 12);
ui.conn_dist_value = uicontrol('Style', 'Text', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.92 0.90 0.06 0.03], ...
    'String', ' ', 'FontSize', 10);
ui.conn_dist = uicontrol('Style', 'Slider', 'Parent', ui.connview, 'Units', 'Normalized', 'Position', [0.91 0.58 0.08 0.32], ...
    'Callback', @conn_dist_Callback);

% %% circlemap
% ui.circlemap = figure('Visible', 'Off', 'Name', 'Surface View', 'NumberTitle', 'Off', 'Position', ...
%     [viewpos(1)-600,viewpos(2),600,600], 'CloseRequestFcn', @view_CloseRequestFcn);
%% surface view window
ui.surfview = figure('Visible', 'Off', 'Name', 'Surface View', 'NumberTitle', 'Off', 'Position', ...
    [viewpos(1)-600,viewpos(2),600,600], 'CloseRequestFcn', @view_CloseRequestFcn);
ui.surfaxes = axes('Parent', ui.surfview, 'NextPlot', 'Add');
axis(ui.surfaxes,'off');

ui.stimelec1_text = uicontrol('Parent', ui.surfview,'Units','normalized','Style','text',...
    'String',' ','FontSize',12,'Position',[0,0.96,0.36,0.03]);

ui.stimelec2_text = uicontrol('Parent', ui.surfview,'Units','normalized','Style','text',...
    'String',' ','FontSize',12,'Position',[0,0.93,0.36,0.03]);

ui.AP = uicontrol('Parent', ui.surfview,'Units','normalized','Style','pushbutton','String','A>P',...
    'Position',[0.01,0.03,0.08,0.05],...
    'Tag','AP',...
    'Callback',@switch_view_Callback);

ui.PA = uicontrol('Parent', ui.surfview,'Units','normalized','Style','pushbutton','String','P>A',...
    'Position',[0.10,0.03,0.08,0.05],...
    'Tag','PA',...
    'Callback',@switch_view_Callback);

ui.SI = uicontrol('Parent', ui.surfview,'Units','normalized','Style','pushbutton','String','S>I',...
    'Position',[0.19,0.03,0.08,0.05],...
    'Tag','SI',...
    'Callback',@switch_view_Callback);
ui.IS = uicontrol('Parent', ui.surfview,'Units','normalized','Style','pushbutton','String','I>S',...
    'Position',[0.28,0.03,0.08,0.05],...
    'Tag','IS',...
    'Callback',@switch_view_Callback);
ui.LR = uicontrol('Parent', ui.surfview,'Units','normalized','Style','pushbutton','String','L>R',...
    'Position',[0.37,0.03,0.08,0.05],...
    'Tag','LR',...
    'Callback',@switch_view_Callback);
ui.RL = uicontrol('Parent', ui.surfview,'Units','normalized','Style','pushbutton','String','R>L',...
    'Position',[0.46,0.03,0.08,0.05],...
    'Tag','RL',...
    'Callback',@switch_view_Callback);

ui.ADT = uicontrol('Parent', ui.surfview,'Units','normalized','Style','pushbutton','String','ADT',...
    'Position',[0.55,0.03,0.12,0.05],...
    'Tag','ADT',...
    'Callback',@ADT_Callback);

ui.opacity_text = uicontrol('Parent', ui.surfview,'Units','normalized','Style','text',...
    'String','Opacity','FontSize',12,'Position',[0.68,0.04,0.08,0.03]);

ui.opacity = uicontrol('Parent', ui.surfview,'Units','normalized','Style','edit',...
    'String','0.5','FontSize',12,...
    'Position',[0.78,0.035,0.05,0.04],'value', 1,...
    'Callback', @opacity_Callback);

ui.activation = uicontrol('Parent', ui.surfview,'Units','normalized','Style','toggle',...
    'String', 'Activation','FontSize', 10, 'Position', [0.84,0.05,0.07,0.03], ...
    'Value', 1, 'Callback', @activation_Callback);

ui.connect = uicontrol('Parent', ui.surfview,'Units','normalized','Style','toggle',...
    'String', 'Connect','FontSize', 10, 'Position', [0.91,0.05,0.07,0.03], ...
    'Value', 0,'Callback', @connect_Callback);

ui.parcellation = uicontrol('Parent', ui.surfview,'Units','normalized','Style','toggle',...
    'String', 'Parcellation','FontSize', 10, 'Position', [0.84,0.02,0.07,0.03], ...
    'Value', 0,'Callback', @parcellation_Callback);

ui.reset = uicontrol('Parent', ui.surfview,'Units','normalized','Style','pushbutton',...
    'String', 'Reset','FontSize', 10, 'Position', [0.91,0.02,0.07,0.03], ...
    'Callback', @reset_Callback);





end



%% ui callback funtions

function dirpath_Callback(hObject, eventdata)

global subinfo
mainpath = get(hObject,'String');
subinfo.mainpath = mainpath;

end


function dirbrowse_Callback(varargin)

global ui subinfo pointVOL

if length(varargin)==2
    mainpath = uigetdir('','Please choose the subject directory');
elseif length(varargin)==1
    mainpath=varargin{1};
end

if mainpath ~=0
    subinfo.mainpath = mainpath;
    set(ui.dirpath,'String',[mainpath filesep])
    
    
    fprintf('Loading relevant files... \n')
    
    fprintf('Loading subject''s infomation...\n')
    if exist([mainpath filesep 'subjectinfo.txt'],'file')
        fid = fopen([mainpath filesep 'subjectinfo.txt']);
        tempinfo = textscan(fid,'%s');
        fclose(fid);
        tempinfo = tempinfo{1};
        subinfo.name = tempinfo{2};
        set(ui.name,'String',subinfo.name)
        subinfo.sex = tempinfo{4};
        if strcmp(subinfo.sex,'male')
            set(ui.sex,'Value',2)
        else set(ui.sex,'Value',3)
        end
        subinfo.age = tempinfo{6};
        set(ui.age,'String',num2str(subinfo.age))
        subinfo.sphere = tempinfo{8};
        if strcmp(subinfo.sphere,'lh')
            set(ui.sphere,'Value',2)
        else if strcmp(subinfo.sphere,'rh')
                set(ui.sphere,'Value',3)
            else set(ui.sphere,'Value',4)
            end
        end
        subinfo.num_elec = str2num(tempinfo{10});
        subinfo.chan_array = cellfun(@str2num,tempinfo(12:12+str2num(tempinfo{10})-1));
        subinfo.chan_array = subinfo.chan_array';
        subinfo.cum_chan_array = cumsum(subinfo.chan_array);
        subinfo.numchan = sum(subinfo.chan_array);
        set(ui.numchan,'String',num2str(subinfo.num_elec));
        set(ui.chan_array,'String',num2str(subinfo.chan_array));
        temp_name = cell2mat(tempinfo(12+str2num(tempinfo{10})+1:end));
        set(ui.chan_ABC,'String',temp_name);
        
        cutinfo = strfind(temp_name, ';');
        for i = 1:length(cutinfo)
            if i == 1
                subinfo.chan_name{i} = temp_name(1:cutinfo(i)-1);
            else
                
                subinfo.chan_name{i} = temp_name(cutinfo(i-1)+1:cutinfo(i)-1);
            end
        end
        
        subinfo.conn_char = cell(subinfo.cum_chan_array(end),1);
        
        n = 0;
        for ien = 1:subinfo.num_elec
            for j = 1:subinfo.chan_array(ien)
                n = n + 1;
                subinfo.conn_char{n} = [subinfo.chan_name{ien} num2str(j)];
            end
        end
    else
        fprintf('Subject Infomation file not found,Please load manually.\n')
    end
    
    %     subinfo.stand_coors = importdata([mainpath filesep 'brain3D' filesep 'MNI152_coordinates_ras.txt']);
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
        subinfo.elec_pairs = elec_pairs;
        
        cum_chan_per_elec = [0 subinfo.cum_chan_array];
        
        for i=1:length(elec_pairs)
            elec_pairs_label{i,1} = [subinfo.conn_char{elec_pairs(i,1)} ,'-', subinfo.conn_char{elec_pairs(i,2)}];
        end
        set(ui.stimelec,'String',elec_pairs_label)
        
    else
        warndlg('Stimulation data missing,please check...')
    end
    %
    fprintf('Loading individual pacellation file... \n')
    
    if exist([subinfo.mainpath filesep 'mri' filesep 'aparc+aseg.nii'],'file')
        mri_info = MRIread([subinfo.mainpath filesep 'mri' filesep 'aparc+aseg.nii']); % input argument
        subinfo.aparc_vol = int32(mri_info.vol); % use int to save space and increase speed
    elseif exist([subinfo.mainpath filesep 'mri' filesep 'aparc+aseg.mgz'],'file')
        mri_info = MRIread([subinfo.mainpath filesep 'mri' filesep 'aparc+aseg.mgz']); % input argument
        subinfo.aparc_vol = int32(mri_info.vol); % use int to save space and increase speed
    else
        fprintf('Failed to load parcellation file,corresponding info may unable to display.\n')
    end
    
    
    if exist([subinfo.mainpath filesep 'brain3D' filesep 'T1.nii'],'file')
        t1info=load_nii([subinfo.mainpath filesep 'brain3D' filesep 'T1.nii']);
        T1vol = t1info.img;
    else
        errordlg('Failed to load T1.nii file.')
    end
    
    % build electrode ROIs
    if exist([subinfo.mainpath filesep 'brain3D' filesep 'coordinates.txt' ],'file') ...
            || exist([subinfo.mainpath filesep 'brain3D' filesep 'autocoordinates.mat' ],'file')
        subinfo.coordinates = ccep_CreateROIs(subinfo);
        subinfo.dist = dist(subinfo.coordinates');
        distance = subinfo.dist;
        subinfo.dist = subinfo.dist - diag(diag(subinfo.dist));
    end
    if  ~exist([subinfo.mainpath filesep 'brain3D' filesep 'distance.mat' ],'file')
        save([subinfo.mainpath filesep 'brain3D' filesep 'distance.mat' ], 'distance');
    end
    % load electrode ROIs
    if exist([subinfo.mainpath filesep 'brain3D' filesep 'elec_image.nii'],'file')
        elecImginfo=load_nii([subinfo.mainpath filesep 'brain3D' filesep 'elec_image.nii']);
        T1vol = T1vol + elecImginfo.img;
    end
    
    if ~isfield(t1info.hdr.hist,'old_affine')
        t1info.hdr.hist.old_affine =  [t1info.hdr.hist.srow_x;...
            t1info.hdr.hist.srow_y; t1info.hdr.hist.srow_z;0 0 0 1];
    end
    
    subinfo.vox2ras = t1info.hdr.hist.old_affine;
    subinfo.ras2vox = inv(subinfo.vox2ras);
    subinfo.vox2ras_tkr = [-1,0,0,128;0,0,1,-128;0,-1,0,128;0,0,0,1];
    subinfo.ras2vox_tkr = inv(subinfo.vox2ras_tkr);
    
    
    ui.Image = single(T1vol);
    ui.size = size(ui.Image);
    pointVOL = round([ui.size(1),ui.size(2),ui.size(3)]/2);
    
    ui.minvalue = min(ui.Image(:));
    ui.min_con = ui.minvalue;
    ui.maxvalue = max(ui.Image(:));
    ui.max_con = ui.maxvalue;
    
    % Creat orthogonal views and cursor marker
    axes(ui.mainaxisc);
    
    ui.Cview = imshow(rot90(squeeze(ui.Image(:,pointVOL(2),:))),[ui.minvalue 0.8*ui.maxvalue]); % Coronal View
    hold on
    ui.vMarkerC = line([0 0],[0 0],'Color',[0 0 1]);
    ui.hMarkerC = line([0 0],[0 0],'Color',[0 0 1]);
    ui.stiC=scatter([],[],'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1]);
    axis(ui.mainaxisc,'equal');
    axis(ui.mainaxisc,'off');
    xlim(ui.mainaxisc,[0 ui.size(1)+1])
    ylim(ui.mainaxisc,[0 ui.size(3)+1])
    
    axes(ui.mainaxiss);
    ui.Sview = imshow(rot90(squeeze(ui.Image(pointVOL(1),:,:))),[ui.minvalue 0.8*ui.maxvalue]); % Sagittal View
    hold on
    ui.vMarkerS = line([0 0],[0 0],'Color',[0 0 1]);
    ui.hMarkerS = line([0 0],[0 0],'Color',[0 0 1]);
    ui.stiS=scatter([],[],'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1]);
    axis(ui.mainaxiss,'equal');
    axis(ui.mainaxiss,'off');
    xlim(ui.mainaxiss,[0 ui.size(2)+1])
    ylim(ui.mainaxiss,[0 ui.size(3)+1])
    
    axes(ui.mainaxisa);
    ui.Aview = imshow(rot90(squeeze(ui.Image(:,:,pointVOL(3)))),[ui.minvalue 0.8*ui.maxvalue]); % AxialCoronal View
    hold on
    ui.vMarkerA = line([0 0],[0 0],'Color',[0 0 1]);
    ui.hMarkerA = line([0 0],[0 0],'Color',[0 0 1]);
    ui.stiA=scatter([],[],'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1]);
    axis(ui.mainaxisa,'equal');
    axis(ui.mainaxisa,'off');
    xlim(ui.mainaxisa,[0 ui.size(1)+1])
    ylim(ui.mainaxisa,[0 ui.size(2)+1])
    
    
    
    
    fprintf('Complete. \n')
end
end



function name_Callback(hObject, eventdata)

global subinfo
str = get(hObject,'String');
subinfo.name = str;

end


function sex_Callback(hObject, eventdata)

global subinfo
val = get(hObject,'Value');
if val == 1
    subinfo.sex = 'male';
else subinfo.sex = 'female';
end

end


function age_Callback(hObject, eventdata)

global subinfo
str = get(hObject,'String');
subinfo.age = str2num(str);

end


function sphere_Callback(hObject, eventdata)

global subinfo
val = get(hObject,'Value');

if val == 1
    subinfo.sphere = 'lh';
else if val == 2
        subinfo.sphere = 'rh';
    else
        subinfo.sphere = 'bi';
    end
end

end


function numchan_Callback(hObject, eventdata)

global subinfo
str = get(hObject,'String');
subinfo.numchan = str2num(str);

end


function chan_array_Callback(hObject, eventdata)

global subinfo
str = get(hObject,'String');
subinfo.chan_array = str2num(str);
subinfo.cum_chan_array = cumsum(subinfo.chan_array);
subinfo.num_elec = length(subinfo.chan_array);
subinfo.numchan = sum(subinfo.chan_array);


set(ui.numchan,'String',num2str(subinfo.numchan));


end

function chan_ABC_Callback(hObject, eventdata)

global subinfo ui

set(ui.chan_ABC,'String',num2str(subinfo.numchan));


end

function stimelec_Callback(hObject, eventdata)

global ui subinfo pointVOL
feval(@stimelec3D_Callback,hObject,eventdata) %surview
val = get(ui.stimelec,'Value');
str = get(ui.stimelec,'String');
% Set current data to the selected data set.

stimelec = subinfo.elec_pairs(val,:);
elec_label = str{val,1};
sepmarker = find(elec_label == '-');
label1 = elec_label(1:sepmarker-1);
label2 = elec_label(sepmarker+1:end);
% [region,percent] = find_roi(subinfo.vol_coordinate(stimelec,:),subinfo.aparc_vol,ui.aparc_aseg,ui.aparc_indx,subinfo.roi_radius);
% region1 = deblank(region{1});
% region2 = deblank(region{1});
% set(ui.stimelec1_text, 'String', [char(450) label1 ':' region1 '  ' num2str(round(100*percent(1))) '%'])
% set(ui.stimelec2_text, 'String', [char(450) label2 ':' region2 '  ' num2str(round(100*percent(2))) '%'])

% plot ccep per a pair electrodes

if exist([subinfo.mainpath filesep 'stimulationdata' filesep 'ccep_elec_' label1 '_' label2 '_All.mat'],'file')
    load([subinfo.mainpath filesep 'stimulationdata' filesep 'ccep_elec_' label1 '_' label2 '_All.mat']);
else
    errordlg('Failed to load corresponding stimulation file.\n');
end

peak=subinfo.connmat(strcmp(subinfo.conn_char,label1),:);
label1=find(strcmp(subinfo.conn_char,label1));
label2=find(strcmp(subinfo.conn_char,label2));

[~,subinfo.stiVOL]=ccep_mapping_final(subinfo,peak(:),[label1,label2]);


%    load and display ccep map
if exist([subinfo.mainpath filesep 'brain3D' filesep 'sccep_map_raw.nii'],'file')
    mapinfo=load_nii([subinfo.mainpath filesep 'brain3D' filesep 'sccep_map_raw.nii']);
    mapvol = mapinfo.img;
else
    errordlg('Failed to load sccep_map_raw.nii file.')
end

ui.ccepImage = single(mapvol);
ui.ccepColormap = jet(round(max(ui.ccepImage(:))));%get colorbar by this matrix


Cccep = rot90(squeeze(ui.ccepImage(:,round(pointVOL(2)),:)));
Cccep_rgb = zeros(ui.size(3),ui.size(1),3);
for i = 1:ui.size(3)
    for j = 1:ui.size(1)
        if Cccep(i,j)~=0
            Cccep_rgb(i,j,:) = ui.ccepColormap(round(Cccep(i,j)),:);
        else
            Cccep_rgb(i,j,:) = [0.4 0.4 0.4];
        end
    end
end

Sccep = rot90(squeeze(ui.ccepImage(round(pointVOL(1)),:,:)));
Sccep_rgb = zeros(ui.size(3),ui.size(2),3);
for i = 1:ui.size(3)
    for j = 1:ui.size(2)
        if  Sccep(i,j)~=0
            Sccep_rgb(i,j,:) = ui.ccepColormap(round(Sccep(i,j)),:);
        else
            Sccep_rgb(i,j,:) = [0.4 0.4 0.4];
        end
        
    end
end

Accep = rot90(squeeze(ui.ccepImage(:,:,ui.size(3)+1-round(pointVOL(3)))));
Accep_rgb = zeros(ui.size(2),ui.size(1),3);
for i = 1:ui.size(2)
    for j = 1:ui.size(1)
        if  Accep(i,j)~=0
            Accep_rgb(i,j,:) = ui.ccepColormap(round(Accep(i,j)),:);
        else
            Accep_rgb(i,j,:) = [0.4 0.4 0.4];
        end
        
    end
end
if isfield(ui,'Cccep')
    set(ui.Cccep,'CData',0)
end
if isfield(ui,'Sccep')
    set(ui.Sccep,'CData',0);
end
if isfield(ui,'Accep')
    set(ui.Accep,'CData',0);
end



perc = 0.4;
axes(ui.mainaxisc);hold on;
ui.Cccep = imshow(Cccep_rgb); % Parcellation Coronal View
alpha(ui.Cccep,perc)

axes(ui.mainaxiss);hold on;
ui.Sccep = imshow(Sccep_rgb); % Parcellation Sagittal View
alpha(ui.Sccep,perc)

axes(ui.mainaxisa);hold on;
ui.Accep = imshow(Accep_rgb); % Parcellation Axial View
alpha(ui.Accep,perc)


end

function stimelec3D_Callback(hObject, eventdata)

global ui subinfo
val = get(ui.stimelec,'Value');
str = get(ui.stimelec,'String');
% Set current data to the selected data set.

stimelec = subinfo.elec_pairs(val,:);
elec_label = str{val,1};
sepmarker = find(elec_label == '-');
label1 = elec_label(1:sepmarker-1);
label2 = elec_label(sepmarker+1:end);

% Show the percentage for belonging to each region
% [region,percent] = find_roi(subinfo.vol_coordinate(stimelec,:),subinfo.aparc_vol,ui.aparc_aseg,ui.aparc_indx,subinfo.roi_radius);
% region1 = deblank(region{1});
% region2 = deblank(region{1});
% set(ui.stimelec1_text, 'String', [char(450) label1 ':' region1 '  ' num2str(round(100*percent{1})) '%'])
% set(ui.stimelec2_text, 'String', [char(450) label2 ':' region2 '  ' num2str(round(100*percent{2})) '%'])


% plot ccep per a pair electrodes

if exist([subinfo.mainpath filesep 'stimulationdata' filesep 'ccep_elec_' label1 '_' label2 '_All.mat'],'file')
    load([subinfo.mainpath filesep 'stimulationdata' filesep 'ccep_elec_' label1 '_' label2 '_All.mat']);
    
    latency_N1_ori = latency_N1;
    latency_P1_ori = latency_P1;
else
    errordlg('Failed to load corresponding stimulation file.\n');
end


%%
%respnse value
peak=subinfo.connmat(strcmp(subinfo.conn_char,label1),:);

% set connectivity matrix cross
try
    %     x = stimelec(1);
    %     y = stimelec(2);
    %     set(ui.conn_cross,'Xdata',x,'Ydata',y,'Color',[0 1 0])
    %     connectivity = subinfo.connmat(y,x);
    %     xlabel = subinfo.conn_char(x);
    %     ylabel = subinfo.conn_char(y);
    %     xind = num2str(x+1-find(subinfo.conn_char(x) == subinfo.conn_char, 1 ));
    %     yind = num2str(y+1-find(subinfo.conn_char(y) == subinfo.conn_char, 1 ));
    %
    %     set(ui.conn_index_x,'String',[xlabel xind])
    %     set(ui.conn_index_y,'String',[ylabel yind])
    %     set(ui.conn_value,'String',num2str(connectivity))
catch
end
%% draw the stimulate electrode pair
try
    %     set(ui.stim_elec_selectc,'XData',subinfo.tkr_coordinate(stimelec,1),'YData',subinfo.tkr_coordinate(stimelec,3), ...
    %         'SizeData', ui.spoint, 'MarkerFaceColor', 'g', 'linewidth', 1);
    %     set(ui.stim_elec_selects,'XData',subinfo.tkr_coordinate(stimelec,2),'YData',subinfo.tkr_coordinate(stimelec,3), ...
    %         'SizeData', ui.spoint, 'MarkerFaceColor', 'g', 'linewidth', 1);
    %     set(ui.stim_elec_selecta,'XData',subinfo.tkr_coordinate(stimelec,1),'YData',subinfo.tkr_coordinate(stimelec,2), ...
    %         'SizeData', ui.spoint, 'MarkerFaceColor', 'g', 'linewidth', 1);
    set(ui.stim_elec_select,'XData',subinfo.tkr_coordinate(stimelec,1),'YData',subinfo.tkr_coordinate(stimelec,2),'ZData',subinfo.tkr_coordinate(stimelec,3), ...
        'SizeData', ui.spoint, 'MarkerEdgeColor', [0 0 0], 'linewidth', 3);
    %     set(ui.stim_elec_text1,'Position',subinfo.tkr_coordinate(stimelec(1),:)+0.5,'String',char(450))
    %     set(ui.stim_elec_text2,'Position',subinfo.tkr_coordinate(stimelec(2),:)+0.5,'String',char(450))
catch
end


%% only the specific indicator

subinfo.peak = peak;


if isfield(ui,'ADT_bar')
    subinfo.first = [0 subinfo.cum_chan_array(1:end-1)]+1;
    subinfo.firstpeak = subinfo.peak;
    subinfo.firstpeak(setdiff(1:length(subinfo.peak),subinfo.first)) = 0;
    subinfo.firstpeak(intersect(1:length(subinfo.peak),subinfo.first)) = 0.8;
    subinfo.otherpeak = subinfo.peak;
    subinfo.otherpeak(subinfo.first) = 0;
    set(ui.ADT_bar(1),'XData',1:length(subinfo.otherpeak))
    set(ui.ADT_bar(1),'YData',subinfo.otherpeak);
    set(ui.ADT_bar(2),'XData',1:length(subinfo.firstpeak))
    set(ui.ADT_bar(2),'YData',subinfo.firstpeak);
    
end

hold on
if isfield(ui,'stim_elec')
    delete(ui.stim_elec);
end
if isfield(ui,'surfview_axes')
    delete(ui.surfview_axes);
end

artifact_elecs=find(peak==0);
val_elecs = setdiff(1:sum(subinfo.chan_array),unique([stimelec artifact_elecs]));

peak(peak==0)=nan;
cvalue= jet(101);
value=min(peak):(max(peak)-min(peak))/100.:max(peak);


[sortedPeak, indxPeak]= sort(peak,'ascend');
PeakIndex = indxPeak(find(sortedPeak>0));


for ielec= val_elecs
    if ismember(ielec,PeakIndex)
        [~,q]=min(abs(peak(ielec)-value));
        ui.stim_elec(ielec,1) = scatter3(subinfo.tkr_coordinate(ielec,1),subinfo.tkr_coordinate(ielec,2),subinfo.tkr_coordinate(ielec,3), ...
            ui.spoint,cvalue(q,:),'fill','Parent',ui.surfaxes,'PickableParts', 'None');
    else
        ui.stim_elec(ielec,1) = scatter3(subinfo.tkr_coordinate(ielec,1),subinfo.tkr_coordinate(ielec,2),subinfo.tkr_coordinate(ielec,3), ...
            ui.spoint,ui.cpoint,'fill','Parent',ui.surfaxes,'PickableParts', 'None');
    end
    
end

ui.surfview_axes = axes('Parent', ui.surfview, 'NextPlot', 'Add','Visible','Off','XTick',[],'YTick',[]);
linkprop([ui.surfaxes,ui.surfview_axes],{'View',...
    'cameraposition',...
    'cameraupvector',...
    'cameratarget',...
    'cameraviewangle'});

if size(cvalue,1)>0 %only set the colorbar
    set(ui.surfview_axes,'Clim',[min(peak) max(peak)]);
    colormap(ui.surfview_axes,cvalue);
end

% ui.stim_colorbar = colorbar(ui.surfview_axes,'Position',[0.8 0.2 0.03 0.6]);
ui.stim_colorbar = colorbar(ui.surfview_axes,'fontsize',15,'fontname','Times New Roman');

set(ui.surfview,'CurrentAxes',ui.surfaxes);
end

function mainwindow_WindowButtonDownFcn(hObject, eventdata)

global ui pointVOL

% set the botton-down function of main window
if strcmp(get(hObject,'SelectionType'),'normal')
    ui.currentaxes = gca;
    set(hObject,'WindowButtonMotionFcn',@wbmcb);
    set(hObject,'WindowButtonUpFcn',@wbucb);
    switch ui.currentaxes
        case ui.mainaxisc
            
            pt = get(gca,'CurrentPoint');
            x = pt(1,1)-0.5; % mouse click coordinate
            y = pt(1,2)-2.5;
            
            if x > 0.5 && x < ui.size(1)+0.5 && y >0.5 && y <  ui.size(3)+0.5
                pointVOL(1) = x;
                pointVOL(3) = y;
                feval(@updatecursor_Callback,hObject, eventdata);
            end
        case ui.mainaxiss
            
            pt = get(gca,'CurrentPoint');
            y = pt(1,1)-0.5; % mouse click coordinates modified
            z = pt(1,2)-2.5;
            if y > 0.5 && y < ui.size(2)+0.5 && z >0.5 && z <  ui.size(3)+0.5
                pointVOL(2) = y;
                pointVOL(3) = z;
                feval(@updatecursor_Callback,hObject, eventdata);
            end
        case ui.mainaxisa
            
            pt = get(gca,'CurrentPoint');
            x = pt(1,1)-0.5; % mouse click coordinates modified
            z = pt(1,2)-2.5;
            if x > 0.5 && x < ui.size(1)+0.5 && z >0.5 && z <  ui.size(2)+0.5
                pointVOL(1) = x;
                pointVOL(2) = ui.size(2)+1-z;
                feval(@updatecursor_Callback,hObject, eventdata);
            end
    end
    
end

    function wbmcb(hObject,eventdata)
        switch ui.currentaxes
            case ui.mainaxisc
                
                pt = get(gca,'CurrentPoint');
                x = pt(1,1)-0.5; % mouse click coordinate
                y = pt(1,2)-2.5;
                
                if x > 0.5 && x < ui.size(1)+0.5 && y >0.5 && y <  ui.size(3)+0.5
                    pointVOL(1) = x;
                    pointVOL(3) = y;
                    feval(@updatecursor_Callback,hObject, eventdata);
                end
            case ui.mainaxiss
                
                pt = get(gca,'CurrentPoint');
                y = pt(1,1)-0.5; % mouse click coordinates modified
                z = pt(1,2)-2.5;
                if y > 0.5 && y < ui.size(2)+0.5 && z >0.5 && z <  ui.size(3)+0.5
                    pointVOL(2) = y;
                    pointVOL(3) = z;
                    feval(@updatecursor_Callback,hObject, eventdata);
                end
            case ui.mainaxisa
                
                pt = get(gca,'CurrentPoint');
                x = pt(1,1)-0.5; % mouse click coordinates modified
                z = pt(1,2)-2.5;
                if x > 0.5 && x < ui.size(1)+0.5 && z >0.5 && z <  ui.size(2)+0.5
                    pointVOL(1) = x;
                    pointVOL(2) = ui.size(2)+1-z;
                    feval(@updatecursor_Callback,hObject, eventdata);
                end
        end
    end

    function wbucb(hObject,eventdata)
        set(hObject,'WindowButtonMotionFcn','');
        set(hObject,'WindowButtonUpFcn','');
    end

end



function mainwindow_WindowKeyPressFcn(hObject, eventdata, handles)

global ui pointVOL

% set the key press function of main window
if isfield(ui,'Image')
    switch gca
        case ui.mainaxisc
            key = get(hObject,'CurrentKey');
            switch key
                case 'leftarrow'
                    if pointVOL(1) > 1.5
                        pointVOL(1) = pointVOL(1)-1;
                    end
                case 'rightarrow'
                    if pointVOL(1) < ui.size(1)-0.5
                        pointVOL(1) = pointVOL(1)+1;
                    end
                case 'uparrow'
                    if pointVOL(3) > 1.5
                        pointVOL(3) = pointVOL(3)-1;
                    end
                case 'downarrow'
                    if pointVOL(3) < ui.size(3)-0.5
                        pointVOL(3) = pointVOL(3)+1;
                    end
                case 'w'
                    if pointVOL(2) < ui.size(2)-0.5
                        pointVOL(2) = pointVOL(2)+1;
                    end
                case 's'
                    if pointVOL(2) > 1.5
                        pointVOL(2) = pointVOL(2)-1;
                    end
            end
            feval(@updatecursor_Callback,hObject, eventdata);
        case ui.mainaxiss
            key = get(hObject,'CurrentKey');
            switch key
                case 'leftarrow'
                    if pointVOL(2) > 1.5
                        pointVOL(2) = pointVOL(2)-1;
                    end
                case 'rightarrow'
                    if pointVOL(2) < ui.size(2)-0.5
                        pointVOL(2) = pointVOL(2)+1;
                    end
                case 'uparrow'
                    if pointVOL(3) > 1.5
                        pointVOL(3) = pointVOL(3)-1;
                    end
                case 'downarrow'
                    if pointVOL(3) < ui.size(3)-0.5
                        pointVOL(3) = pointVOL(3)+1;
                    end
                case 'w'
                    if pointVOL(1) < ui.size(1)-0.5
                        pointVOL(1) = pointVOL(1)+1;
                    end
                case 's'
                    if pointVOL(1) > 1.5
                        pointVOL(1) = pointVOL(1)-1;
                    end
            end
            feval(@updatecursor_Callback,hObject, eventdata);
        case  ui.mainaxisa
            key = get(hObject,'CurrentKey');
            switch key
                case 'leftarrow'
                    if pointVOL(1) > 1.5
                        pointVOL(1) = pointVOL(1)-1;
                    end
                case 'rightarrow'
                    if pointVOL(1) < ui.size(1)-0.5
                        pointVOL(1) = pointVOL(1)+1;
                    end
                case 'uparrow'
                    if pointVOL(2) < ui.size(2)-0.5
                        pointVOL(2) = pointVOL(2)+1;
                    end
                case 'downarrow'
                    if pointVOL(2) > 1.5
                        pointVOL(2) = pointVOL(2)-1;
                    end
                case 'w'
                    if pointVOL(3) < ui.size(3)-0.5
                        pointVOL(3) = pointVOL(3)-1;
                    end
                case 's'
                    if pointVOL(3) > 1.5
                        pointVOL(3) = pointVOL(3)+1;
                    end
            end
            feval(@updatecursor_Callback,hObject, eventdata);
    end
    
end

end


function mainwindow_CloseRequestFcn(hObject,eventdata)

global ui subinfo
delete(hObject);


if isfield(ui,'connview')
    try
        delete(ui.connview)
    catch
    end
end

if isfield(ui,'surfview')
    try
        delete(ui.surfview)
    catch
    end
end

if isfield(ui,'circleview')
    try
        delete(ui.circleview)
    catch
    end
end

ui = [];
subinfo = [];
clear ui subinfo
clc
fprintf('Thanks for using MRIESviewer. \n')
end



% --- Executes on cursor reposition.
function updatecursor_Callback(hObject, eventdata)
% update current cursor location and the corresponding coordinates
global ui subinfo pointVOL

if isfield(ui,'Image')
    pointvol1 = [ui.size(1)+1-round(pointVOL(1)),round(pointVOL(3)),round(pointVOL(2))];
    if isfield(subinfo,'stiVOL')
        
        stiVOL = mean (subinfo.stiVOL);
    else
        stiVOL=[0 0 0];
    end
    if abs(stiVOL(3)-pointvol1(3))<4
        ui.xC=ui.size(1)+1-round(stiVOL(1));
        ui.yC=stiVOL(2);
    else
        ui.xC=[];
        ui.yC=[];
    end
    
    if abs(stiVOL(1)-pointvol1(1))<4
        ui.xS=stiVOL(3);
        ui.yS=stiVOL(2);
    else
        ui.xS=[];
        ui.yS=[];
    end
    
    
    if abs(stiVOL(2)-pointvol1(2))<4
        ui.xA=ui.size(1)+1-round(stiVOL(1));
        ui.yA=ui.size(2)+1-stiVOL(3);
    else
        ui.xA=[];
        ui.yA=[];
    end
    
    pointRAS = subinfo.vox2ras*[pointvol1-1 1]';
    pointRAS = pointRAS(1:3)';
    pointRAS = round(pointRAS*100)/100;
    
    set(ui.vMarkerC, 'XData', [pointVOL(1) pointVOL(1)], 'YData', [pointVOL(3)-ui.size(3)-1 pointVOL(3)+ui.size(3)+1]);
    set(ui.hMarkerC, 'XData', [pointVOL(1)-ui.size(1)-1 pointVOL(1)+ui.size(1)+1], 'YData', [pointVOL(3) pointVOL(3)]);
    set(ui.vMarkerS, 'XData', [pointVOL(2) pointVOL(2)], 'YData', [pointVOL(3)-ui.size(3)-1 pointVOL(3)+ui.size(3)+1]);
    set(ui.hMarkerS, 'XData', [pointVOL(2)-ui.size(2)-1 pointVOL(3)+ui.size(3)+1], 'YData', [pointVOL(3) pointVOL(3)]);
    set(ui.vMarkerA, 'XData', [pointVOL(1) pointVOL(1)], 'YData', [-pointVOL(2) 2*ui.size(2)+2-pointVOL(2)]);
    set(ui.hMarkerA, 'XData', [pointVOL(1)-ui.size(1)-1 pointVOL(1)+ui.size(1)+1], 'YData', [ui.size(2)+1-pointVOL(2) ui.size(2)+1-pointVOL(2)]);
    
    set(ui.Cview,'CData',rot90(squeeze(ui.Image(:,round(pointVOL(2)),:)))); % Set Coronal View
    set(ui.Sview,'CData',rot90(squeeze(ui.Image(round(pointVOL(1)),:,:)))); % Set Sagittal View
    set(ui.Aview,'CData',rot90(squeeze(ui.Image(:,:,ui.size(3)+1-round(pointVOL(3)))))); % Set Axial View
    
    
    
    if isfield(ui,'Cccep')
        Cccep = rot90(squeeze(ui.ccepImage(:,round(pointVOL(2)),:)));
        Cccep_rgb = zeros(ui.size(3),ui.size(1),3);
        for i = 1:ui.size(3)
            for j = 1:ui.size(1)
                if Cccep(i,j)~=0
                    Cccep_rgb(i,j,:) = ui.ccepColormap(round(Cccep(i,j)),:);
                else
                    Cccep_rgb(i,j,:) = [0.4 0.4 0.4];
                end
            end
        end
        
        Sccep = rot90(squeeze(ui.ccepImage(round(pointVOL(1)),:,:)));
        Sccep_rgb = zeros(ui.size(3),ui.size(2),3);
        for i = 1:ui.size(3)
            for j = 1:ui.size(2)
                if Sccep(i,j)~=0
                    Sccep_rgb(i,j,:) = ui.ccepColormap(round(Sccep(i,j)),:);
                else
                    Sccep_rgb(i,j,:) = [0.4 0.4 0.4];
                end
            end
        end
        
        Accep = rot90(squeeze(ui.ccepImage(:,:,ui.size(3)+1-round(pointVOL(3)))));
        Accep_rgb = zeros(ui.size(2),ui.size(1),3);
        for i = 1:ui.size(2)
            for j = 1:ui.size(1)
                if Accep(i,j)~=0
                    Accep_rgb(i,j,:) = ui.ccepColormap(round(Accep(i,j)),:);
                else
                    Accep_rgb(i,j,:) = [0.4 0.4 0.4];
                end
            end
        end
        
        set(ui.Cccep,'CData',Cccep_rgb)
        set(ui.Sccep,'CData',Sccep_rgb)
        set(ui.Accep,'CData',Accep_rgb)
        set(ui.overlay_value,'String',num2str(ui.ccepImage(round(pointVOL(1)),round(pointVOL(2)),ui.size(3)+1-round(pointVOL(3)))))
        
    end
    
    set(ui.stiC,'XData',ui.xC,'YData',ui.yC);
    uistack(ui.stiC,'top')
    set(ui.stiS,'XData',ui.xS,'YData',ui.yS);
    uistack(ui.stiS,'top')
    set(ui.stiA,'XData',ui.xA,'YData',ui.yA);
    uistack(ui.stiA,'top')
    set(ui.volcoord,'String',num2str(pointvol1))
    set(ui.rascoord,'String',num2str(pointRAS))
    set(ui.underlay_value,'String',num2str(ui.Image(round(pointVOL(1)),round(pointVOL(2)),ui.size(3)+1-round(pointVOL(3)))))
    
    
    
end

end



function circle_button_Callback(hObject,eventdata)

global ui subinfo


if isempty(subinfo.connmat)
    if exist([subinfo.mainpath filesep 'Matrix' filesep 'conn_matrix_' subinfo.method '_.mat'],'file')
        fprintf('Loading connectivity matrix... \n')
        a = load([subinfo.mainpath filesep 'Matrix' filesep 'conn_matrix_' subinfo.method '_.mat']);
        c = fieldnames(a);
        subinfo.connmat = getfield(a,c{1});
    elseif exist([subinfo.mainpath filesep 'stimulationdata'],'dir')
        fprintf('Calculating connectivity matrix.. \n')
        subinfo.connmat = calculate_connectivity_matrix(subinfo);
    else
        warndlg('Loading stimulation data failed,please chech the corresponding files.')
    end
end

if isempty(subinfo.connthresh)
    subinfo.connthresh=0;
end

R1=10;
R2=9;

[roicortex,percent] = find_roi(subinfo.vol_coordinate,subinfo.aparc_vol,ui.aparc_aseg,ui.aparc_indx,subinfo.roi_radius);
roicortex = deblank(roicortex);
[cortex_unique,Iall,Iuni] = unique(roicortex);
cortex_new = cortex_unique;

for n = 1:length(cortex_unique)
    cortex_elecnum(n) = length(find(Iuni==n));
    cortex_elec{n} = find(Iuni == n);
end
m = 0;
for n = 1:length(cortex_unique)
    if ~ismember(cortex_unique(n),ui.ccep_cortex)
        cortex_new(n-m) = [];
        cortex_elecnum(n-m) = [];
        cortex_elec(n-m) = [];
        m = m+1;
    end
end
textstr = char(cortex_new);
textstr(find(textstr=='_')) = '-';

cvalue = cool(length(cortex_unique));
interval = 2*pi/(sum(cortex_elecnum)+length(cortex_new));
delta = interval;
for i = 1:length(cortex_new)
    deltap(i,1) = delta;
    for j = 1:cortex_elecnum(i)
        alpha = linspace(0,interval,100)+delta;
        x1 = R1*cos(alpha);
        y1 = R1*sin(alpha);
        x2 = R2*cos(alpha);
        y2 = R2*sin(alpha);
        x = [x2 fliplr(x1)];
        y = [y2 fliplr(y1)];
        plot(x,y,'Parent',ui.circleaxes1);
        hold on
        ui.circle_elec = fill(x,y,cvalue(ceil(length(cortex_unique)/2)),'Parent',ui.circleaxes1,'Tag', num2str(cortex_elec{i}(j)),'ButtonDownFcn', @circle_elec_ButtonDownFcn);
        rotation_angle = mean(alpha)*180/pi;
        if rotation_angle>90 && rotation_angle<270
            ui.circle_elec_text(cortex_elec{i}(j)) = text((R1+0.2)*cos(mean(alpha)),(R1+0.2)*sin(mean(alpha)),' ', 'rotation',rotation_angle-180,...
                'HorizontalAlignment','right','Parent',ui.circleaxes2, 'Visible', 'Off', 'PickableParts', 'None');
        else
            ui.circle_elec_text(cortex_elec{i}(j)) = text((R1+0.2)*cos(mean(alpha)),(R1+0.2)*sin(mean(alpha)),' ', 'rotation',rotation_angle,...
                'HorizontalAlignment','left','Parent',ui.circleaxes2, 'Visible', 'Off', 'PickableParts', 'None');
        end
        delta = max(alpha);
    end
    deltap(i,2) = delta;
    delta = delta+interval;
end
%     axis square;
for i = 1:length(cortex_new)
    position(i,:) = [R2*cos(mean(deltap(i,:))),R2*sin(mean(deltap(i,:)))];
    textxp = (R1+0.5)*cos(mean(deltap(i,:)));
    textyp = (R1+0.5)*sin(mean(deltap(i,:)));
    rotation_angle = mean(deltap(i,:))*180/pi;
    if rotation_angle>90&&rotation_angle<270
        ui.circle_text(i) =text(textxp,textyp,textstr(i,:),'HorizontalAlignment','right','rotation',rotation_angle-180,'fontsize',20,'Parent',ui.circleaxes1);
    else
        ui.circle_text(i) = text(textxp,textyp,textstr(i,:),'HorizontalAlignment','left','rotation',rotation_angle,'fontsize',20,'Parent',ui.circleaxes1);
    end
end

direction = zeros(length(cortex_new));
for stim_cortex = 1:length(cortex_new)
    stim_elec = cortex_elec{stim_cortex};
    for resp_cortex = 1:length(cortex_new)
        if resp_cortex ~= stim_cortex
            resp_elec = cortex_elec{resp_cortex};
            direction_mark = 0;
            conn_temp = subinfo.connmat(stim_elec,resp_elec);
            conn_temp(conn_temp<=subinfo.connthresh) = 0;
            if sum(find(conn_temp>0))>=1
                if direction(resp_cortex,stim_cortex) == 1;
                    direction_mark = 1;  %% bidirectional
                end
                ui.circle_map= bulge_byx(position(stim_cortex,1),position(stim_cortex,2), ...
                    position(resp_cortex,1),position(resp_cortex,2),2,direction_mark,ui.circleaxes2);
                direction(stim_cortex,resp_cortex) = 1;
            end
        end
    end
end
colormap(ui.circleaxes2,jet(32))
colormap(ui.circleaxes1,gray(32))
if strcmpi(get(hObject,'String'),'circle map')
    if strcmpi(get(ui.circleview,'Visible'),'On')
        set(ui.circleview,'Visible','Off')
    else
        set(ui.circleview,'Visible','On')
    end
end



end

function circle_elec_ButtonDownFcn(hObject,eventdata)

global ui subinfo
elec_index = str2double(get(hObject,'Tag'));
label = subinfo.conn_char(elec_index);
% ind = num2str(elec_index+1-find(subinfo.conn_char(elec_index) == subinfo.conn_char, 1 ));
ind = num2str(elec_index);
elec_label = strcat(label,num2str(ind));
if strcmp(get(ui.circle_elec_text(elec_index),'String'),' ')
    [region,percent] = find_roi(subinfo.vol_coordinate(elec_index,:),subinfo.aparc_vol,ui.aparc_aseg,ui.aparc_indx,subinfo.roi_radius);
    region = deblank(region{1});
    set(ui.circle_elec_text(elec_index), 'String', [elec_label ': ' num2str(round(100*percent{1})) '%'],'Color',[1 0 0],'fontsize',14);
end

if strcmpi(get(ui.circle_elec_text(elec_index),'Visible'),'off')
    set(ui.circle_elec_text(elec_index),'Visible','on')
else
    set(ui.circle_elec_text(elec_index),'Visible','off')
end

end

function coordinate_Callback(hObject,eventdata)

global ui subinfo
fprintf(['loading contacts loaction file...'])
raw_coors = [];
[filename, pathname, ~] = uigetfile('*.mat;*.txt','Please select the electrode coordinates file');
if filename ~= 0
    if strcmp(filename(end-2:end),'txt')
        subinfo.raw_coors = importdata([pathname filename]);
    else
        if strcmp(filename(end-2:end),'mat')
            a = load([pathname filename]);
            c = fieldnames(a);
            subinfo.raw_coors = getfield(a,c{1});
            
        end
    end
    
    if isempty(subinfo.chan_array) && size(subinfo.raw_coors,2) == 5
        table = tabulate(floor(subinfo.raw_coors(:,2)/100));
        subinfo.chan_array = table(:,2)';
        subinfo.cum_chan_array = cumsum(subinfo.chan_array);
        subinfo.num_elec = length(subinfo.chan_array);
        subinfo.numchan = sum(subinfo.chan_array);
        
        set(ui.numchan,'String',num2str(subinfo.numchan));
        set(ui.chan_array,'String',num2str(subinfo.chan_array));
    end
    
    ras_coordinate = subinfo.raw_coors(:,(end-2:end));
    subinfo.tkr_coordinate = subinfo.vox2ras_tkr*subinfo.ras2vox*[ras_coordinate ones(size(ras_coordinate,1),1)]';
    subinfo.tkr_coordinate = subinfo.tkr_coordinate(1:3,:)';
    subinfo.vol_coordinate = subinfo.ras2vox*[ras_coordinate ones(size(ras_coordinate,1),1)]';
    subinfo.vol_coordinate = subinfo.vol_coordinate(1:3,:)';
    
    if subinfo.num_elec<9
        elec_label = char((1:subinfo.num_elec)'+64);
    elseif subinfo.num_elec>=9
        elec_label = char([1:8 10:subinfo.num_elec+1]'+64);
    end
    all_pointxyz = mat2cell(subinfo.tkr_coordinate, subinfo.chan_array, 3);
    
    % point size and color
    ui.spoint = 30;
    ui.cpoint = [ 0 0 0];
    

    if ~isfield(ui,'array')
        for ielec = 1:size(subinfo.tkr_coordinate,1)
            ui.array(ielec) = scatter3(subinfo.tkr_coordinate(ielec,1),subinfo.tkr_coordinate(ielec,2),subinfo.tkr_coordinate(ielec,3), 30, [1 1 1],'filled',...
                'Parent',ui.surfaxes, 'Tag', num2str(ielec), 'ButtonDownFcn', @surf_elec_ButtonDownFcn);
            if any(ielec==subinfo.cum_chan_array)
                ind=find(ielec==subinfo.cum_chan_array);
                
                ui.array_text(ind) = text(subinfo.tkr_coordinate(ielec,1)+5,subinfo.tkr_coordinate(ielec,2)+5,subinfo.tkr_coordinate(ielec,3)+5,subinfo.chan_name{ind}, ...
                    'fontsize',30,'FontWeight','demi' , 'Parent',ui.surfaxes, 'Visible', 'on', 'PickableParts', 'None');
            end
        end
    else
        for ielec = 1:size(subinfo.tkr_coordinate,1)
            set(ui.array(ielec), 'XData', subinfo.tkr_coordinate(ielec,1), 'YData', subinfo.tkr_coordinate(ielec,2), 'ZData', subinfo.tkr_coordinate(ielec,3));
            if any(ielec==subinfo.cum_chan_array)
                ind=find(ielec==subinfo.cum_chan_array);
                set(ui.array_text(ind), 'Position', subinfo.tkr_coordinate(ielec,:)+5, 'PickableParts', 'None');
                
            end
        end
    end
    
    subinfo.dist = squareform(pdist(subinfo.tkr_coordinate(:,1:3),'euclidean'));
    
end
fprintf(['success.\n'])
end


function connbutton_Callback(hObject,eventdata)

global ui subinfo

if strcmpi(get(ui.connview,'Visible'),'On')
    set(ui.connview,'Visible','Off')
else
    set(ui.connview,'Visible','On')

end

if ~isfield(ui,'connmatrix')
    if isempty(subinfo.connmat)
        if exist([subinfo.mainpath filesep 'Matrix' filesep 'conn_matrix_' subinfo.method '.mat'],'file')
            fprintf('Loading connectivity matrix... \n')
            a = load([subinfo.mainpath filesep 'Matrix' filesep 'conn_matrix_' subinfo.method '.mat']);
            c = fieldnames(a);
            subinfo.connmat = getfield(a,c{1});
        elseif exist([subinfo.mainpath filesep 'stimulationdata'],'dir')
            fprintf('Calculating connectivity matrix... \n')
            subinfo.connmat = calculate_connectivity_matrix(subinfo);
        else
            warndlg('Loading stimulation data failed,please chech the corresponding files. \n')
        end
    end
    ui.connmatrix = imagesc(subinfo.connmat, 'Parent', ui.connaxes,'ButtonDownFcn', @connmat_ButtonDownFcn);

    
    
    chan_name_j = cell(subinfo.num_elec,1);
    for ien = 1:subinfo.num_elec
        %         indxj = strfind(subinfo.chan_name{ien},';');
        %         chan_name_j{ien} = subinfo.chan_name{ien}(1:indxj-1);
        chan_name_j{ien} = subinfo.chan_name{ien};
    end
    
    set(ui.connaxes,'ytick',[1 subinfo.cum_chan_array(1:end-1)+1],'yticklabel',chan_name_j,'tickdir','out','YDir','normal')
    set(ui.connaxes,'xtick',[1 subinfo.cum_chan_array(1:end-1)+1],'xticklabel',chan_name_j,'tickdir','out')
    %     text(-3,-3,'A','Parent',ui.connaxes)
    xlabel(ui.connaxes, 'Record');ylabel(ui.connaxes, 'Stim')
    grid(ui.connaxes, 'On')
    hold(ui.connaxes)
    ui.conn_cross = plot(0,0,'+','Color',[1 1 1],'LineWidth',1,'MarkerSize',8,'Parent',ui.connaxes);
    subinfo.distthresh = 0;
    subinfo.connthresh = 0;
end


end


function connmat_ButtonDownFcn(hObject,eventdata)

global ui subinfo

pt = get(gca,'CurrentPoint');
x = round(pt(1,1));
y = round(pt(1,2));
set(ui.conn_cross,'Xdata',pt(1,1)-0.1,'Ydata',pt(1,2)+0.4,'Color',[0 1 0])
connectivity = subinfo.connmat(y,x);
xlabel = subinfo.conn_char(x);
ylabel = subinfo.conn_char(y);

set(ui.conn_index_r,'String',xlabel{1})
set(ui.conn_index_s,'String',ylabel{1})
set(ui.conn_value,'String',num2str(connectivity))

str = get(ui.stimelec,'String');
% stim_index = strncmp(str,[ylabel yind],length([ylabel yind]));
stim_index = strncmp(str,ylabel(1),length(ylabel{1}));
val = find(stim_index == 1);
set(ui.stimelec,'Value',val)
% stimelec_Callback(hObject,eventdata);

end

function method_Callback(hObject, eventdata)

global ui subinfo
val = get(hObject,'Value');
str = get(hObject,'String');
subinfo.method = str{val};

if exist([subinfo.mainpath filesep 'Matrix' filesep 'conn_matrix_' subinfo.method '.mat'],'file')
    fprintf(['Loading ', str{val} ,' connectivity  matrix... '])
    a = load([subinfo.mainpath filesep 'Matrix' filesep 'conn_matrix_' subinfo.method '.mat']);
    c = fieldnames(a);
    subinfo.connmat = getfield(a,c{1});
else
    subinfo.connmat = calculate_connectivity_matrix(subinfo);
end
if isfield(ui,'connmatrix')
    connectivity = subinfo.connmat;
    set(ui.connmatrix,'CData',connectivity);
end

if strcmpi(get(ui.surfview,'Visible'),'on')
    if isfield(ui,'stim_elec')
        delete(ui.stim_elec);
        ui = rmfield(ui,'stim_elec');
        feval(@stimelec_Callback,hObject,eventdata)
    end
end


if strcmpi(get(ui.circleview,'Visible'),'on')
    
    if isfield(ui,'circle_map')
        delete(ui.circle_map)
        delete(ui.circle_text)
        ui = rmfield(ui,{'circle_map','circle_text'});
        feval(@circle_button_Callback,hObject,eventdata)
    else
        feval(@circle_button_Callback,hObject,eventdata)
    end
    
end

fprintf(['success.\n'])
end

function conn_thresh_Callback(hObject,eventdata)

global ui subinfo
percent = get(hObject,'Value');
connectivity = subinfo.connmat;
subinfo.connthresh = percent*max(subinfo.connmat(:));
connectivity(connectivity < subinfo.connthresh) = 0;
set(ui.conn_thresh_value,'String',num2str(subinfo.connthresh))

distmat = subinfo.dist;
distmat(distmat < subinfo.distthresh) = 0;
distmat(distmat>0) = 1;
subinfo.conn_dist = connectivity.*distmat;
set(ui.connmatrix,'CData',subinfo.conn_dist)



if strcmpi(get(ui.circleview,'Visible'),'on')
    
    if isfield(ui,'circle_map')
        delete(ui.circle_map)
        delete(ui.circle_text)
        delete(ui.circle_colorbar);
        ui = rmfield(ui,{'circle_map','circle_text','circle_colorbar'});
        feval(@circle_button_Callback,hObject,eventdata)
    else
        feval(@circle_button_Callback,hObject,eventdata)
    end
    
end


end

function conn_dist_Callback(hObject,eventdata)

global ui subinfo
connectivity = subinfo.connmat;
connectivity(connectivity < subinfo.connthresh) = 0;
percent = get(hObject,'Value');
distmat = subinfo.dist;

subinfo.distthresh = percent*max(distmat(:));
distmat(distmat < subinfo.distthresh) = 0;
distmat(distmat>0) = 1;
set(ui.conn_dist_value,'String',num2str(round(subinfo.distthresh)))

subinfo.conn_dist = connectivity.*distmat;
set(ui.connmatrix,'CData',subinfo.conn_dist)
end


function surfacebutton_Callback(hObject,eventdata)

global ui subinfo

if strcmpi(get(ui.surfview,'Visible'),'On')
    set(ui.surfview,'Visible','Off')
else set(ui.surfview,'Visible','On')
    
    if ~isfield(ui,'surf')
        try
            if strcmp(subinfo.sphere,'lh') || strcmp(subinfo.sphere,'rh')
                [vertices, faces]=read_surf([subinfo.mainpath filesep 'surf' filesep subinfo.sphere '.pial']);
                subinfo.faces = faces+1;
                subinfo.vertices = vertices;
                ui.surf = patch(struct(...
                    'vertices', subinfo.vertices, 'faces', subinfo.faces), ...
                    'Parent',ui.surfaxes, ...
                    'FaceColor',[230,228,216]./255, ...
                    'FaceAlpha',0.5, ...
                    'EdgeColor', 'none', ...
                    'PickableParts', 'None');
                axis(ui.surfaxes,'equal')
                axis(ui.surfaxes,'off')
                ui.light = camlight;
                set(ui.light,'position',[0 0 2000])
                lighting(ui.surfaxes,'gouraud')
                material dull
                rotate3d on
            else
                [vertex_coords1, faces1]=read_surf([subinfo.mainpath filesep 'surf' filesep 'lh.pial']);
                faces1 = faces1 + 1;
                [vertex_coords2, faces2]=read_surf([subinfo.mainpath filesep 'surf' filesep 'rh.pial']);
                faces2 = faces2 + 1;
                subinfo.stereonum = [size(vertex_coords1,1) size(faces1,1)];
                subinfo.faces = [faces1;faces2+size(vertex_coords1,1)];
                subinfo.vertices = [vertex_coords1;vertex_coords2];
                ui.surf = patch(struct(...
                    'vertices', subinfo.vertices, 'faces', subinfo.faces), ...
                    'Parent',ui.surfaxes, ...
                    'FaceColor',[230,228,216]./255, ...
                    'FaceAlpha',0.5, ...
                    'EdgeColor', 'none', ...
                    'PickableParts', 'None');
                axis(ui.surfaxes,'equal')
                axis(ui.surfaxes,'off')
                ui.light = camlight;
                set(ui.light,'position',[0 0 2000])
                lighting(ui.surfaxes,'gouraud')
                material dull
                rotate3d on
                
            end
            
            ui.stim_elec_select = scatter3(0,0,0,1,'Parent',ui.surfaxes,'MarkerEdgeColor','w','linewidth',1);
            ui.stim_elec_text1 = text(0,0,0,' ','Parent',ui.surfaxes,'Color',[0 1 0],'FontSize',15,'PickableParts', 'None');
            ui.stim_elec_text2 = text(0,0,1,' ','Parent',ui.surfaxes,'Color',[0 1 0],'FontSize',15,'PickableParts', 'None');
            %             ui.stim_colorbar = colorbar(ui.surfaxes,'TickLength',0,'Position',[0.9 0.4 0.03 0.5]);
            
        catch
            fprintf('Undefined subject infomation,surface file not found. \n')
        end
    end
    
end

end


function surf_elec_ButtonDownFcn(hObject,eventdata)

global ui subinfo
elec_index = str2double(get(hObject,'Tag'));

if strcmp(get(ui.array_text(elec_index),'String'),' ')
    [region,percent] = find_roi(subinfo.vol_coordinate(elec_index,:),subinfo.aparc_vol,ui.aparc_aseg,ui.aparc_indx,subinfo.roi_radius);
    region = deblank(region{1});
    set(ui.array_text(elec_index), 'String', [region '  ' num2str(round(100*percent{1})) '%'],'Color',[1 0 0],'fontsize',14)
end

if strcmpi(get(ui.array_text(elec_index),'Visible'),'off')
    set(ui.array_text(elec_index),'Visible','on')
else
    set(ui.array_text(elec_index),'Visible','off')
end

end


function switch_view_Callback(hObject,eventdata)

global ui
if isfield(ui,'surf')
    switch get(hObject,'Tag')
        case 'AP'
            view(ui.surfaxes,180,0);
            set(ui.light,'position',get(ui.surfaxes,'CameraPosition'))
        case 'PA'
            view(ui.surfaxes,0,0);
            set(ui.light,'position',get(ui.surfaxes,'CameraPosition'))
        case 'SI'
            view(ui.surfaxes,0,90);
            set(ui.light,'position',get(ui.surfaxes,'CameraPosition'))
        case 'IS'
            view(ui.surfaxes,0,-90);
            set(ui.light,'position',get(ui.surfaxes,'CameraPosition'))
        case 'LR'
            view(ui.surfaxes,-90,0);
            set(ui.light,'position',get(ui.surfaxes,'CameraPosition'))
        case 'RL'
            view(ui.surfaxes,90,0);
            set(ui.light,'position',get(ui.surfaxes,'CameraPosition'))
    end
end

end


function ADT_Callback(hObject,eventdata)

global ui subinfo
if ~isfield(ui,'ADT_axis')
    ui.ADT_axis = axes('Parent', ui.surfview,'Units','normalized', 'Position',[0.7,0.092,0.28,0.24], ...
        'XTick', [], 'XTickLabel', [], ...
        'ButtonDownFcn', @ADT_axis_ButtonDownFcn);
    
elseif strcmpi(get(ui.ADT_axis,'Visible'),'off')
    try
        set(ui.ADT_axis,'Visible','on')
        set(ui.ADT_bar,'Visible','on')
        set(ui.ADT_line,'Visible','on')
    end
    
elseif strcmpi(get(ui.ADT_axis,'Visible'),'on')
    try
        set(ui.ADT_axis,'Visible','off')
        set(ui.ADT_bar,'Visible','off')
        set(ui.ADT_line,'Visible','off')
    end
    
end

if ~isfield(ui,'ADT_bar') && ~isempty(subinfo.peak)
    subinfo.first = [0 subinfo.cum_chan_array(1:end-1)]+1;
    subinfo.firstpeak = subinfo.peak;
    subinfo.firstpeak(setdiff(1:length(subinfo.peak),subinfo.first)) = 0;
    subinfo.firstpeak(intersect(1:length(subinfo.peak),subinfo.first)) = 0.8;
    subinfo.otherpeak = subinfo.peak;
    subinfo.otherpeak(subinfo.first) = 0;
    ui.ADT_bar(1) = bar(ui.ADT_axis,subinfo.otherpeak,'FaceColor','b','EdgeColor','b');hold on;
    ui.ADT_bar(2) = bar(ui.ADT_axis,subinfo.firstpeak,'FaceColor','r','EdgeColor','r');
    ui.ADT_line = plot(ui.ADT_axis,[0 256],[0 0]);
    set(ui.ADT_axis,'XLim',[1 length(subinfo.peak)])
    set(ui.ADT_axis,'XTick', [], 'XTickLabel', [])
    set(ui.ADT_axis,'ButtonDownFcn', @ADT_axis_ButtonDownFcn)
    
elseif isfield(ui,'ADT_bar') && ~isempty(subinfo.peak)
    subinfo.first = [0 subinfo.cum_chan_array(1:end-1)]+1;
    subinfo.firstpeak = subinfo.peak;
    subinfo.firstpeak(setdiff(1:length(subinfo.peak),subinfo.first)) = 0;
    subinfo.firstpeak(intersect(1:length(subinfo.peak),subinfo.first)) = 0.8;
    subinfo.otherpeak = subinfo.peak;
    subinfo.otherpeak(subinfo.first) = 0;
    set(ui.ADT_bar(1),'XData',1:length(subinfo.otherpeak))
    set(ui.ADT_bar(1),'YData',subinfo.otherpeak);
    set(ui.ADT_bar(2),'XData',1:length(subinfo.firstpeak))
    set(ui.ADT_bar(2),'YData',subinfo.firstpeak);
    set(ui.ADT_axis,'XLim',[1 length(subinfo.peak)])
    set(ui.ADT_axis,'XTick', [], 'XTickLabel', [])
    set(ui.ADT_axis,'ButtonDownFcn', @ADT_axis_ButtonDownFcn)
    
end

end


function ADT_axis_ButtonDownFcn(hObject,eventdata)

global ui subinfo

pt = get(gca,'CurrentPoint');
subinfo.thresh_percent = pt(1,2)+0.02;
activation = get(ui.activation,'Value');
if subinfo.thresh_percent >= 0 && subinfo.thresh_percent <= 1 && activation == 1
    
    try
        set(ui.ADT_line,'YData',[subinfo.thresh_percent subinfo.thresh_percent])
    end
    if isfield(ui,'stim_elec')
        for i = 1:size(ui.stim_elec,1)
            try
                activation = get(ui.stim_elec(i,1),'SizeData');
                if activation < 400*subinfo.thresh_percent
                    set(ui.stim_elec(i,1),'Visible','Off')
                    
                else  set(ui.stim_elec(i,1),'Visible','On')
                end
                
            end
        end
    end
    
end
end


function opacity_Callback(hObject,eventdata)

global ui
if isfield(ui,'surf')
    opacity = str2num(get(hObject,'String'));
    if opacity >= 0 && opacity <= 1
        set(ui.surf,'FaceAlpha',opacity)
    else
        warndlg('Please input a transparency value between 0 and 1.')
        
    end
end

end


function activation_Callback(hObject,eventdata)

global ui
stat = get(hObject,'Value');

if isfield(ui,'stim_elec')
    switch stat
        case 0
            for i = 1:size(ui.stim_elec,1)
                try
                    set(ui.stim_elec(i,1),'Visible','Off')
                catch
                end
                
            end
        case 1
            for i = 1:size(ui.stim_elec,1)
                try
                    set(ui.stim_elec(i,1),'Visible','On')
                catch
                end
            end
            
    end
    
end

end

function connect_Callback(hObject,eventdata)

global ui subinfo

if isfield(ui,'array') && isfield(subinfo,'connmat')
    
    if isfield(ui,'connect_line')
        delete(ui.connect_line)
    end
    
    if get(hObject,'Value') == 1
        
        connectivity = subinfo.connmat;
        % Display connectivity lines at least half maximum of connectivity
        if subinfo.connthresh < max(subinfo.connmat(:))/2
            connthresh = max(subinfo.connmat(:))/2;
        else
            connthresh = max(subinfo.connmat(:))/4;
        end
        connectivity(connectivity < connthresh) = 0;
        [m,k]=find(connectivity);
        for iline = 1:length(m)
            ui.connect_line(iline,:) = plot_3Dcurve_colordir(subinfo.tkr_coordinate(m(iline),:),subinfo.tkr_coordinate(k(iline),:),2,[ones(11,1) zeros(11,2)]);
        end
        
        
    end
    
    
end

end


function parcellation_Callback(hObject,eventdata)

global ui subinfo
stat = get(hObject,'Value');

if isfield(ui,'surf')
    switch stat
        case 0
            try
                set(ui.surf,'FaceColor',[230,228,216]./255)
            catch
            end
        case 1
            
            if strcmp(subinfo.sphere,'lh') || strcmp(subinfo.sphere,'rh')
                
                [~, L, ct] = read_annotation([subinfo.mainpath filesep 'label' filesep subinfo.sphere '.aparc.annot']);
                
            else
                [~, L1, ct1] = read_annotation([subinfo.mainpath filesep 'label' filesep 'lh.aparc.annot']);
                [~, L2, ct2] = read_annotation([subinfo.mainpath filesep 'label' filesep 'rh.aparc.annot']);
                L = [L1;L2];
                ct.numEntries = ct1.numEntries+ct2.numEntries;
                ct.table = cat(1,ct1.table,ct2.table);
            end
            
            cdata  = repmat([0.5 0.5 0.5], size(subinfo.faces,1),1);
            
            for iroi = 2:ct.numEntries
                roi_indx = iroi;
                roi_vertex_indx = find(L==ct.table(roi_indx,5));
                
                % find the faces including at least one of three vertices
                roi_face_indx=find(sum(ismember(subinfo.faces,roi_vertex_indx),2));
                
                % update the cdata for all rois
                cdata(roi_face_indx,:) = repmat(ct.table(roi_indx,1:3)./255,[size(roi_face_indx,1) 1]);
                set(ui.surf,'facecolor','flat','facevertexcdata',cdata)
            end
            
    end
    
end

end


function reset_Callback(hObject,eventdata)

global ui subinfo

zoom(ui.surfview,'out')

set(ui.activation,'Value',1)
feval(@activation_Callback,ui.activation,eventdata)

set(ui.parcellation,'Value',0)
feval(@parcellation_Callback,ui.parcellation,eventdata)

set(ui.opacity,'String','0.5')
feval(@opacity_Callback,ui.opacity,eventdata)

set(ui.connect,'Value',0)
if isfield(ui,'connect_line')
    delete(ui.connect_line)
end


end


function view_CloseRequestFcn(hObject,eventdata)

set(hObject,'Visible','Off')

end


