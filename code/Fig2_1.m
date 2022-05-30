
%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
% 
% Figure 2-1. Individual Target-locked results
%
% Created: Thur 21 Abril 2022, 19:19
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Thur 30 May 2022, 15:32
% Last edited by: Martin Esparza-Iaizzo
% 
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

% Data directory:
DataFolder = '/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/datasets/dataset2_ref2target';

% Change current folder to Data Folder
cd(DataFolder);
fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected subject \n');


% Initialize target variables
subj = 10; 
tpoints = 751; 
total_R_attR = zeros(subj,tpoints);
total_L_attR = zeros(subj,tpoints);
total_R_attL = zeros(subj,tpoints);
total_L_attL = zeros(subj,tpoints);

%% Main loop thru participants

for i = 1:10

% Select subject
str = sprintf('P00%i_S01_ALH_raw.mat',i);
load(str);
data_power_left = data_phase;

str = sprintf('P00%i_S01_ARH_raw.mat',i);
load(str);
data_power_right = data_phase;

fprintf('File loaded! \n');
fprintf('-----------------------------------------------------------------\n');

% Clear vars
    clearvars -except data_power data_power_left data_power_right NumP SubjectInfo i...
        DataAnalysis total_R_attR total_L_attR total_R_attL total_L_attL p_value_mat iter 
  

%% PREPROCESSING OF RAW DATA & TIME-FREQUENCY ANALYSIS

fprintf('Applying preprocessing to raw data... \n');
% Preprocessing

cfg = [];
cfg.detrend     = 'no';        % Removing of linear trends
cfg.demean      = 'yes';       % Setting of the baseline correction
cfg.dftfilter   = 'yes';       % Application of a line noise filter
cfg.dftfreq     = 50;          % Frequency which get excluded as line noise
cfg.hpfilter    = 'yes';       % High-Pass-Filter On
cfg.hpfiltdir   = 'twopass';   % Forward filtering.
cfg.hpfreq      = 0.16;        % High-Pass Frequency
cfg.hpfilttype  = 'but';       % Filter Type Butterworth (IIR)
cfg.hpfiltord   = 5;           % Filter Order, watch out for Instability!
cfg.lpfilter    = 'yes';       % Low-Pass Filter On
cfg.lpfiltdir   = 'twopass';   % Forward filtering.
cfg.lpfreq      = 45;          % Low-Pass Frequency
cfg.lpfilttype  = 'but';       % Filter Type
cfg.lpfiltord   = 16;          % Filter Order

data_filt_left = ft_preprocessing(cfg, data_power_left);
data_filt_right = ft_preprocessing(cfg, data_power_right);

% Time Frequency Analysis
fprintf('Applying Morlet wavelet anaylisis to preprocessed data... \n');

cfg             = [];
cfg.channel    = {'Fz';'O2';'FC1';'P6';'FC2';'P8';'O1';'P5';'P7'};
cfg.output      = 'fourier';%'powandcsd';
cfg.method      = 'wavelet';
cfg.pad         = 'maxperlen';
cfg.foi         = linspace(9.54, 14.31, 10);
cfg.width       = 5;
cfg.toi         = -0.75:0.002:0.75;
cfg.keeptrials  = 'yes';
cfg.avgoverfreq = 'yes';
power_left      = ft_freqanalysis(cfg, data_filt_left); % Attended Left trials
power_right     = ft_freqanalysis(cfg, data_filt_right); % Attended Right trials

% Connectivity analysis

% Fronto-medial to Parietal left ROI
cfg             = [];
cfg.method      = 'plv';
cfg.channelcmb  = {'F*','O1';'F*','P5';'F*','P7'};
FM_PL_attL = ft_connectivityanalysis(cfg, power_left);
FM_PL_attR = ft_connectivityanalysis(cfg, power_right);

% Fronto-medial to Parietal right ROI
cfg             = [];
cfg.method      = 'plv';
cfg.channelcmb  = {'F*','O2';'F*','P6';'F*','P8'};
FM_PR_attL = ft_connectivityanalysis(cfg, power_left);
FM_PR_attR = ft_connectivityanalysis(cfg, power_right);

% Average over frequencies and channels
FM_PL_attL = squeeze(mean(FM_PL_attL.plvspctrm,1:2));
FM_PL_attR = squeeze(mean(FM_PL_attR.plvspctrm,1:2));
FM_PR_attL = squeeze(mean(FM_PR_attL.plvspctrm,1:2));
FM_PR_attR = squeeze(mean(FM_PR_attR.plvspctrm,1:2));

total_L_attL(i,:) = FM_PL_attL';
total_R_attL(i,:) = FM_PR_attL';
total_L_attR(i,:) = FM_PL_attR';
total_R_attR(i,:) = FM_PR_attR';
end

%% Temporal variables

color_contra = [70 70 70]./255;
color_ipsi = [170 170 170]./255;
shade_ipsi = [220 220 220]./255;

color_attR = [46, 89, 132]./255;
color_attL = [114, 165, 197]./255;
shade_attL = [202, 211, 232]./255;
fontsize = 16;

x1 = [-0.2,-0.1,0.1,0.2]; % Auxiliary variable for x-axis placement

%% Plotting

for subj = 1:9
[~,pre] = find(power_left.time >= -0.2 & power_left.time <= 0);
[~,post] = find(power_left.time >= 0.2 & power_left.time <= 0.4);

violin_pre_ipsi = [total_L_attL(subj,pre) total_R_attR(subj,pre)];
violin_pre_contra = [total_L_attR(subj,pre) total_R_attL(subj,pre)];

violin_post_ipsi = [total_L_attL(subj,post) total_R_attR(subj,post)];
violin_post_contra = [total_L_attR(subj,post) total_R_attL(subj,post)];


f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17.6 20]*1.75);

ax = axes('Position',[0.1 0.7 0.35 0.15]);
ax.PositionConstraint = 'innerposition';

Y{:,1}=violin_pre_ipsi;
Y{:,2}=violin_pre_contra;
Y{:,3}=violin_post_ipsi;
Y{:,4}=violin_post_contra;
scaling_factor = 5.5;

violin(Y,'x',x1.*scaling_factor,'facecolor',...
    [color_ipsi;color_contra;color_ipsi;color_contra],'edgecolor','k',...
    'facealpha',1,'bw',0.04,'mc','k','medc',[]);
xlim([-0.3 0.3].*scaling_factor)
ylim([0 0.7])
yticks([0 0.2 0.4 0.6]); yticklabels({'0',' .2','.4','.6'})
str = sprintf('P0%s – PLV', num2str(subj));
ylabel(str,'fontweight','bold')
xticks([-0.15 0.15].*scaling_factor); xticklabels({'Pre-target','Post-target'})
set(ax,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

%% Save figure ––––––– Uncomment and edit to save to personalised location 
% str = sprintf('s0%s.pdf',num2str(subj));
% exportgraphics(gcf,str,'BackgroundColor','white','ContentType','vector')

end
