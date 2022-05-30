
%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 5-1. Lat. index individual figures. 
%
% Created: Sat 16 Apr 2022, 22:55
% Author: Irene Vigue-Guix
% 
% Last edited:  Mon 30 May 2022, 15:32
% Last edited by: Martin Esparza-Iaizzo
% 
%% Change path to project

% Directory Martin:
ProjectFolder = '/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/analysis';

% Change current folder to Project Folder
cd(ProjectFolder);
fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected project \n');

% PATH SETUP: Set up data analysis pipeline (path folders, toolboxes...)
setupPath_CVSA_dataAnalysis;

%% CHANGE PWD TO SUBJECT'S DATA FOLDER AND LOAD DATASETS

% Add path of target-locked data. 
DataFolder = '/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/datasets/dataset1_ref2cue';
cd(DataFolder);

fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected subject \n');

% Select participant
fprintf('Select a raw data file from a participant: \n');
d = dir;
fn = {d.name};
[nPart,tf] = listdlg('PromptString','Select a file:', 'SelectionMode','single',...
                           'ListString',fn(4:end));  
% Load dataset selected
load(fn{nPart+3}); %load(fn{NumP});
fprintf('File loaded! \n');
fprintf('-----------------------------------------------------------------\n');

% Clear vars
clearvars -except data_power data_power_left data_power_right NumP ...
    DataAnalysis SubjectInfo nsubj nPart

%% LOAD DATA OF INTEREST
% Load subject info
load('SubjectInformation.mat');

% Number Subject
idx_subj = [2 4 8 9 10 13 1 3 5 7];
nsubj = find(idx_subj == nPart);
% Select IAF/IFoI of the participant
IAF = SubjectInfo.IAP_Hz(nPart);
    

%% REFLECTING FILTERED DATA 
% Select data
cfg = [];
cfg.latency = [0 1.5]; % we only have from 0 to 1.5 s
data_left_refl = ft_selectdata(cfg, data_power_left);
data_right_refl = ft_selectdata(cfg, data_power_right);
data_refl = ft_selectdata(cfg, data_power);
% create new time vector for reflecting window
vec_time = {linspace(-1.5,3,2251)};
cell_left_time = repmat(vec_time,length(data_left_refl.time), 1)';
cell_right_time = repmat(vec_time,length(data_right_refl.time), 1)';
cell_time = repmat(vec_time,length(data_refl.time), 1)';

% reflect EEG matrix
cell_left_trial = cell(1,length(data_left_refl.time));
for k = 1:length(data_left_refl.time)
    vec_left_trial = fliplr(data_left_refl.trial{1,k});
    cell_left_trial{k} = [vec_left_trial(:,1:end-1) data_left_refl.trial{1,k} vec_left_trial(:,2:end)];
end
cell_right_trial = cell(1,length(data_right_refl.time));
for k = 1:length(data_right_refl.time)
    vec_right_trial = fliplr(data_right_refl.trial{1,k});
    cell_right_trial{k} = [vec_right_trial(:,1:end-1) data_right_refl.trial{1,k} vec_right_trial(:,2:end)];
end
cell_trial = cell(1,length(data_refl.time));
for k = 1:length(data_refl.time)
    vec_trial = fliplr(data_refl.trial{1,k});
    cell_trial{k} = [vec_trial(:,1:end-1) data_refl.trial{1,k} vec_trial(:,2:end)];
end

% create double-fieldtrip structure
data_left_refl.trial = cell_left_trial;
data_left_refl.time = cell_left_time;
data_right_refl.trial = cell_right_trial;
data_right_refl.time = cell_right_time;
data_refl.trial = cell_trial;
data_refl.time = cell_time;

%% TIME-FREQUENCY ANALYSIS
% Multitaper as in MBC thesis by Irene
fprintf('Applying multitaper analysis to raw data... \n');

%%%
% Time Frequency Analysis
cfg             = [];
cfg.width       = 6;
cfg.output      = 'pow';
cfg.method      = 'wavelet';
cfg.foi         = (IAF-1:1:IAF+1); %IAF; %unique([flip(2.^(log2(IAF):-.25:1)) 2.^(log2(IAF):.25:5.4)]);
cfg.keeptrials  = 'yes';
cfg.pad         = 'nextpow2';
cfg.toi         = -1.5:0.02:3; %0:0.02:1.5;
power_left      = ft_freqanalysis(cfg, data_left_refl);
power_right     = ft_freqanalysis(cfg, data_right_refl);
power           = ft_freqanalysis(cfg, data_refl);



%% DIVISION
% Divide into left/right ROI

% ATTENDED LEFT
% Right ROI
cfg                     = [];
cfg.channel             = {'P6';'P8';'PO4';'O2'}; %{'P6';'P8';'O2'};
cfg.avgoverchan         = 'yes';
cfg.avgoverfreq         = 'yes';
cfg.latency             = [0 1]; %[0 1.5];
data_powerL_rightROI    = ft_selectdata(cfg, power_left);
% Left ROI
cfg                 = [];
cfg.channel         = {'P7';'P5';'PO3';'O1'}; %{'P7';'P3';'O1'}; 
cfg.avgoverchan     = 'yes';
cfg.avgoverfreq     = 'yes';
cfg.latency         = [0 1]; %[0 1.5];
data_powerL_leftROI = ft_selectdata(cfg, power_left);

% ATTENDED RIGHT
% Right ROI
cfg                     = [];
cfg.channel             = {'P6';'P8';'PO4';'O2'}; %{'P4';'P8';'O2'};
cfg.avgoverchan         = 'yes';
cfg.avgoverfreq         = 'yes';
cfg.latency             = [0 1]; %[0 1.5];
data_powerR_rightROI    = ft_selectdata(cfg, power_right);
% Left ROI
cfg                 = [];
cfg.channel         = {'P7';'P5';'PO3';'O1'}; %{'P7';'P3';'O1'};
cfg.avgoverchan     = 'yes';
cfg.avgoverfreq     = 'yes';
cfg.latency         = [0 1]; %[0 1.5];
data_powerR_leftROI = ft_selectdata(cfg, power_right);

%% LATERALIZATION INDEX
% Compute lateralization index based on Thut et al. (2006)

% Attended Right (Att_R)
dp_R_avg = (squeeze(data_powerR_rightROI.powspctrm)+squeeze(data_powerR_leftROI.powspctrm))/2;
dp_R_diff = (squeeze(data_powerR_rightROI.powspctrm)-squeeze(data_powerR_leftROI.powspctrm));
dp_R_LI = (dp_R_diff)./dp_R_avg;
mean_dp_R_LI = mean(dp_R_LI);
std_dp_R_LI = std(dp_R_LI);
meanmean_dp_R_LI = mean(mean_dp_R_LI);
sem_dp_R_LI  = std(dp_R_LI)/sqrt(size(dp_R_LI,2));
timevec = data_powerR_rightROI.time;

% Attended Left (Att_L)
dp_L_avg = (squeeze(data_powerL_rightROI.powspctrm)+squeeze(data_powerL_leftROI.powspctrm))/2;
dp_L_diff = (squeeze(data_powerL_rightROI.powspctrm)-squeeze(data_powerL_leftROI.powspctrm));
dp_L_LI = (dp_L_diff)./dp_L_avg;
mean_dp_L_LI = mean(dp_L_LI);
std_dp_L_LI = std(dp_L_LI);
sem_dp_L_LI  = std(dp_L_LI)/sqrt(size(dp_L_LI,2));
meanmean_dp_L_LI = mean(mean_dp_L_LI);

%% STATS (MEAN LAT INDEX)
% Ttest para mean values of lat index for att. left & right
[h,p,ci,stats] = ttest2(mean_dp_L_LI(~isnan(mean_dp_L_LI))',...
    mean_dp_R_LI(~isnan(mean_dp_R_LI))', 'Alpha', 0.05, 'Tail', 'left');
% Effect size (Cohen's d)
d = Cohen_d(mean_dp_L_LI(~isnan(mean_dp_L_LI))', ...
    mean_dp_R_LI(~isnan(mean_dp_R_LI))','unpaired');

%% STATS (LAT INDEX OVER TIME)
% (5) UNPAIRED T-TEST one-tailed (left - X less than Y)
C1 = dp_L_LI; 
C2 = dp_R_LI;
[h,p,ci,stats] = ttest2(C1, C2, 'Alpha', 0.05, 'Tail','left'); 
alpha=0.05; % set alpha level

% EFFECT SIZE (Cohen's d)
d = Cohen_d(mean(C1),mean(C2),'unpaired');
    
%% Plotting

color_attR = [46, 89, 132]./255;
color_attL = [114, 165, 197]./255;
shade_attL = [202, 211, 232]./255;
fontsize = 18;
scaling_factor = 7.5;
x1 = [0, 0.1];

%%%%%% --------

f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17.6 20]*1.75);

ax1 = axes('Position',[0.1 0.725 0.225 0.175]);
ax1.PositionConstraint = 'innerposition';
clear Y
Y{:,1} = mean_dp_L_LI(~isnan(mean_dp_L_LI));
Y{:,2} = mean_dp_R_LI(~isnan(mean_dp_R_LI));
violin(Y,'x',x1.*scaling_factor,'facecolor',...
    [color_attL;color_attR],'edgecolor','k',...
    'facealpha',1,'bw',0.09,'mc','k','medc',[]); hold on;
% scatter(0.3.*(rand(1,10) - 0.5) + ones(1,10).*x1(1).*scaling_factor,mean_ind_values_AttL,50,...
%     'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.15); hold on;
% scatter(0.3.*(rand(1,10) - 0.5) + ones(1,10).*x1(2).*scaling_factor,mean_ind_values_AttR,50,...
%     'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.15);
ylim([-0.4 0.8]);
xlim([-0.05 0.15].*scaling_factor);
% yticks([-0.4 0 0.4 0.8]); yticklabels({'-.4','0','.4','.8'})
yticks([]); yticklabels({''});
xticks([0 0.1].*scaling_factor); xticklabels({'Att. Left','Att. Right'})
str = sprintf('P%s – Lat. Index                    ', num2str(nsubj));
ylabel(str,'fontweight','bold')
set(ax1,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

ax2 = axes('Position',[0.1 0.62 0.225 0.075]);
ax2.PositionConstraint = 'innerposition';
stdshade(dp_R_LI,0.25,'-',color_attR,color_attR,timevec+0.5); hold on;
stdshade(dp_L_LI,0.5,'-',color_attL,shade_attL,timevec+0.5); hold on;
ylim([-0.2 0.6])
if nsubj == 1 % P01 (1 DE LOS 10)
    v=axis; % get current axis limits
    plot(timevec+0.5,(p<=alpha)*100-100+v(3)*(.95),'k.','MarkerSize',7.5); % PLOT P-VALUES
end
xticks([0.5 1 1.5]); xticklabels({'0.5','1.0','1.5'})
% yticks([0 0.5]); yticklabels({'0','.5'})
yticks([]); yticklabels({''});
xlabel('Time (s)')
set(ax2,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')
% exportgraphics(gcf,'irene1.pdf','BackgroundColor','white','ContentType','vector')
marker_size = 10;

%% Save figure ––––––– Uncomment and edit to save to personalised location 
% str = sprintf('s%s.pdf',num2str(nsubj));
% cd('/Users/martinesparzaiaizzo/Desktop/@martirene/Figures/Figs_temp/Extended data/Figure 5-1')
% % exportgraphics(gcf,str,'BackgroundColor','white','ContentType','vector')
