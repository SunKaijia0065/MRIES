% Convert data from edf format to mat as the edf_match.txt
% updated by Kaijia Sun

function ccep_edf2mat_ALL(varargin)
if isempty(varargin)
    datapath = uigetdir('','Please choose the subject directory');
else
    datapath=varargin{1};
end

edfPath=[datapath,filesep,'edf',filesep];

files = dir([edfPath filesep  '*.edf']);

fileNames={files.name};
edfMPath=[datapath,filesep,'edf',filesep,'edf_match.txt'];
fid = fopen(edfMPath);
formatSpec = '%s %s  %s';
dict = textscan(fid,formatSpec,'delimiter','\t');           
fclose(fid);
    for i =1:length(fileNames)
        filename=fileNames{i};
        [~, nameEdf, ext] = fileparts(filename);
        index=find(ismember(dict{1},nameEdf));
        edfName=filename;
        dictEle=dict{2}(index);
        dictCon=dict{3}(index);
        ccep_edf2mat_bh(edfName,edfPath,dictEle{1},dictCon{1})


    end



end

