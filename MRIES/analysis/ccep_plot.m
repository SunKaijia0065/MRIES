function ccep_plot
%
% plot ccep curve of the specified stimulated electrodes
%
close all

datapath = uigetdir('','Please choose the Patient stimulationdata directory');
files = dir([datapath filesep 'ccep*_All.mat']);

% [FileName,PathName,~] = uigetfile('*.mat','Please choose the interested stimulated electrodes');
cd(datapath)
for ifile = 1:length(files)
    filename = files(ifile).name;
    load([datapath filesep filename])
    indx = strfind(filename,'_');
    
    elec1 = filename(indx(end-2)+1:indx(end-1)-1);
    elec2 = filename(indx(end-1)+1:indx(end)-1);
    
    load('../subj_elec_info.mat')
    
    ScSz = get( groot, 'Screensize');
    
    pos_grid =  [ScSz(1)+30 ScSz(2)+30  ScSz(3)-30  ScSz(4)-30];
    grid1 = figure('Name',['ccep-elec-' elec1 '-' elec2]);
    set(grid1,'Position', ScSz,'color',[1 1 1],'MenuBar','figure');
    
    nrow = length(numelec);
    ncol = max(numelec);
    col_ind = 1:ncol:nrow*ncol;
    e_trans = [];
    for ielec = 1:nrow
        e_trans = [e_trans col_ind(ielec)+numelec(ielec)-1:-1:col_ind(ielec)];
    end
    
    gridPos = DivideScreen(nrow,ncol,ScSz,40,60,pos_grid);
    val_elecs=setdiff(total_elecs,trigelec);
    for ielec = val_elecs
        axes('Position',gridPos{e_trans(ielec)});
        if ismember(ielec,stimelec)
            fprintf(['elec #' num2str(ielec) ': STIM ' chan_name{ielec} '\n']  );
        else
            fprintf(['elec #' num2str(ielec) ':' chan_name{ielec} '\n']  );
        end
        if ismember(ielec,stimelec)
            title(['STIM: ' chan_name{ielec} ] ,'color','r');
            continue
        end
        if ismember(ielec,badelec)
            title (['BAD:' chan_name{ielec}],'color','b');
            continue
        end
        
        plot(-win(1)+1/Fs:1/Fs:win(2),ccep1(:,ielec),'c');
        hold on
        plot(-win(1)+1/Fs:1/Fs:win(2),ccep2(:,ielec),'k');
        ccep = [ccep1(:,ielec) ccep2(:,ielec)];
        if abs(min(min(ccep))) > max(max(ccep))
            ylim([min(min(ccep))-10 -min(min(ccep))+10]);
        else
            ylim([-max(max(ccep))-10 max(max(ccep))+10]);
        end
        
        %     if abs(min(min(ccep))) < 300 && abs(max(max(ccep))) < 300
        %         ylim([-300 300]);
        %     end
        ylim([-1000 1000]);
        ylimit = get(gca,'Ylim');
        %% find the earlier N1 latency and amplitude
        if ~isnan(latency_N1(ielec,1)) && ~isnan(peak_N1(ielec,1))
            if fwhm_dur(ielec,1) < 100
                plot(latency_N1(ielec,1),peak_N1(ielec,1),'rs','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',8)
                
            end
        end
        %% find the earlier P1 latency and amplitude
        if ~isnan(latency_P1(ielec,1)) && ~isnan(peak_P1(ielec,1))
            if fwhm_dur(ielec,2) < 100
                plot(latency_P1(ielec,1),peak_P1(ielec,1),'rs','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8)
            end
            
        end
        
        %%
        plot(-win(1)+1/Fs:1/Fs:win(2),stat.h(:,ielec)*ylimit(2)/2.*sign(stat.z(:,ielec)),'m');
        %     area([40 800], [ylimit(2) ylimit(2)], ylimit(1), 'FaceColor', [0.5 .5 .5]);
        
        
        set(gca,'YDir','reverse')
        xlim([-0.1 0.5]);
        
        if ~any((stat.h(round((win(1)+0.005)*Fs):end,ielec)>0) | (stat.h(round((win(1)+0.005)*Fs):end,ielec)<0))
            %         title (['REJECTED:' chan_name{ielec}],'color','b');
        else
            title (['elec' num2str(ielec) ':' chan_name{ielec}])
        end
        
        box off
        % axcopy
        
        
    end
    hgsave(['ccep_elec_' elec1 '_'  elec2 '.fig']);
    close(grid1)
end

% %% specify the interested electrodes for double check
% val_elecs=input('Please specify the interested electrodes. eg. [103 104]: ');
% % val_elecs = [];
% if ~isempty(val_elecs)
%     for ielec = val_elecs
%         figure('Name',['elec' num2str(ielec) ':' chan_name{ielec}])
%         
%         
%         plot(-win(1)+1/Fs:1/Fs:win(2),ccep1(:,ielec),'c');
%         hold on
%         plot(-win(1)+1/Fs:1/Fs:win(2),ccep2(:,ielec),'k');
%         ccep = [ccep1(:,ielec) ccep2(:,ielec)];
%         if abs(min(min(ccep))) > max(max(ccep))
%             ylim([min(min(ccep))-10 -min(min(ccep))+10]);
%         else
%             ylim([-max(max(ccep))-10 max(max(ccep))+10]);
%         end
%         
%         %     if abs(min(min(ccep))) < 300 && abs(max(max(ccep))) < 300
%         %         ylim([-300 300]);
%         %     end
%         ylim([-500 500]);
%         ylimit = get(gca,'Ylim');
%         %% find the earlier N1 latency and amplitude
%         if ~isnan(latency_N1(ielec,1)) && ~isnan(peak_N1(ielec,1))
%             if fwhm_dur(ielec,1) < 100
%                 plot(latency_N1(ielec,1),peak_N1(ielec,1),'rs','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',8)
%                 
%             end
%         end
%         %% find the earlier P1 latency and amplitude
%         if ~isnan(latency_P1(ielec,1)) && ~isnan(peak_P1(ielec,1))
%             if fwhm_dur(ielec,2) < 100
%                 plot(latency_P1(ielec,1),peak_P1(ielec,1),'rs','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8)
%             end
%             
%         end
%         
%         %%
%         plot(-win(1)+1/Fs:1/Fs:win(2),stat.h(:,ielec)*ylimit(2)/2.*sign(stat.z(:,ielec)),'m');
%         %     area([40 800], [ylimit(2) ylimit(2)], ylimit(1), 'FaceColor', [0.5 .5 .5]);
%         
%         
%         set(gca,'YDir','reverse')
%         xlim([-0.1 0.5]);
%         
%         box off
%         %     axcopy
%         
%         
%     end
% end


