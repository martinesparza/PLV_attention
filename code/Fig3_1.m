
%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 3-1. Cue-locked individual figures. 
%
% Created: Sat 23 Apr 2022, 14:07
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Mon 30 May 2022, 15:32
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

%% Load subjects

% Add cue-locked data
DataFolder_1 = '/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/datasets/dataset1_ref2cue';
% Add raw data
DataFolder_0 = '/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/datasets/dataset0_rawdata';
addpath(DataFolder_0);
addpath(DataFolder_1);

cd(DataFolder_0)

% Initialize target variables
subj = 10; 
tpoints = 151; 
total_R_attR = zeros(subj,tpoints);
total_L_attR = zeros(subj,tpoints);
total_R_attL = zeros(subj,tpoints);
total_L_attL = zeros(subj,tpoints);

% Auxiliary mapping variable
subjects = [2 4 13 14 15 18 1 3 5 12];

%% Main loop thru participants

for i = 1:10
    
% Select subject
if subjects(i) < 10
    str = sprintf('OffLine_000%i',subjects(i)); 
    str1 = sprintf('cue_preproc_trials_samelength_offline_000%i.mat',subjects(i));
    load(str1);
end
if subjects(i) >= 10
    str = sprintf('OffLine_00%i',subjects(i));
    str1 = sprintf('cue_preproc_trials_samelength_offline_00%i.mat',subjects(i));
    load(str1);
end
fname=str;
offset=2; % In seconds
cfg=[];

cfg.dataset                 =[fname '.eeg'];
cfg.headerfile              =[fname '.vhdr'];
cfg.eventfile               =[fname '.vmrk'];

hdr   = ft_read_header(cfg.headerfile);
cfg.method='trial';

%Now we get the trial matrix from the data already segmented
clear trl_matrix_left trl_matrix_right

trl_matrix_left(:,1:2)=data_power_left.sampleinfo;
trl_matrix_left(:,3)=offset*hdr.Fs; %This is the offset (where are we putting the zero, in samples)
trl_matrix_left(:,4)=data_power_left.trialinfo(:,1); %This is the kind of trial we have (hit/mis)

trl_matrix_right(:,1:2)=data_power_right.sampleinfo;
trl_matrix_right(:,3)=offset*hdr.Fs; %This is the offset (where are we putting the zero, in samples)
trl_matrix_right(:,4)=data_power_right.trialinfo(:,1); %This is the kind of trial we have (hit/mis)


%% Preprocessing

fprintf('Applying preprocessing to raw data... \n');

% cfg = [];
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

cfg.trl=trl_matrix_left; % Use attemded left trials
data_filt_left = ft_preprocessing(cfg);

for k = 1:length(data_filt_left.time)
    data_filt_left.time{1,k} = data_filt_left.time{1,k} - 2;
end

cfg.trl=trl_matrix_right; % Use attemded left trials
data_filt_right = ft_preprocessing(cfg);

for k = 1:length(data_filt_right.time)
    data_filt_right.time{1,k} = data_filt_right.time{1,k} - 2;
end

%% TIME FREQUENCY ANALYSIS

fprintf('Applying morlet wavelet analysis to raw data... \n');

% Mirrowing

cfg = [];
cfg.latency = [-1.5 3];
data_refl_left = ft_selectdata(cfg, data_filt_left);
data_refl_right = ft_selectdata(cfg, data_filt_right);

% create new time vector for reflecting window
vec_time_left = {linspace(-1.5,3,2251)};
cell_time_left = repmat(vec_time_left,length(data_refl_left.time), 1)';

vec_time_right = {linspace(-1.5,3,2251)};
cell_time_right = repmat(vec_time_right,length(data_refl_right.time), 1)';


% reflect EEG matrix
cell_trial_left = cell(1,length(data_refl_left.time));
for k = 1:length(data_refl_left.time)
    vec_trial = fliplr(data_refl_left.trial{1,k});
    cell_trial_left{k} = [vec_trial(:,1:end-1) data_refl_left.trial{1,k} vec_trial(:,2:end)];
end 

cell_trial_right = cell(1,length(data_refl_right.time));
for k = 1:length(data_refl_right.time)
    vec_trial = fliplr(data_refl_right.trial{1,k});
    cell_trial_right{k} = [vec_trial(:,1:end-1) data_refl_right.trial{1,k} vec_trial(:,2:end)];
end 

% create double-fieldtrip structure
data_refl_left.trial = cell_trial_left;
data_refl_left.time = cell_time_left;

data_refl_right.trial = cell_trial_right;
data_refl_right.time = cell_time_right;

% Time-frequency analysis
cfg             = [];
cfg.output      = 'fourier';
cfg.channel    = {'Fz';'O2';'FC1';'P6';'FC2';'P8';'O1';'P5';'P7'};
cfg.method      = 'wavelet';
cfg.pad         = 'maxperlen';
cfg.foi         = linspace(9.54, 14.31, 5);
cfg.width       = 5;
cfg.toi         = -0.25:0.01:1.25; % Relevant time is from 0 to 1s
cfg.keeptrials  = 'yes';
power_left      = ft_freqanalysis(cfg, data_refl_left);
power_right     = ft_freqanalysis(cfg, data_refl_right);

% Time windows

[~,first] = find(power_left.time >= 0 & power_left.time < 0.2);
[~,second] = find(power_left.time >= 0.2 & power_left.time < 0.4);
[~,third] = find(power_left.time >= 0.4 & power_left.time < 0.6);
[~,fourth] = find(power_left.time >= 0.6 & power_left.time < 0.8);
[~,fifth] = find(power_left.time >= 0.8 & power_left.time < 1);

% Fronto-medial to Parietal left ROI
cfg             = [];
cfg.method      = 'plv';
cfg.channelcmb  = {'F*','O1';'F*','P5';'F*','P7';'F*','PO3'};
FM_PL_attL = ft_connectivityanalysis(cfg, power_left);
FM_PL_attR = ft_connectivityanalysis(cfg, power_right);

% Fronto-medial to Parietal right ROI
cfg             = [];
cfg.method      = 'plv';
cfg.channelcmb  = {'F*','O2';'F*','P6';'F*','P8';'F*','PO4'};
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

%% Plotting

color_contra = [70 70 70]./255;
color_ipsi = [170 170 170]./255;
shade_ipsi = [220 220 220]./255;

color_attL = [204, 153, 204]./255;
color_attR = [102, 153, 204]./255;
fontsize = 16;
scaling_factor = 7.5;

% Auxiliary plotting variable
x1 = [0.07, 0.13, 0.27, 0.33, 0.47, 0.53, 0.67, 0.73, 0.87, 0.93];

cc_ipsi(:,:,1) = total_R_attR;
cc_ipsi(:,:,2) = total_L_attL;

cc_contra(:,:,1) = total_R_attL;
cc_contra(:,:,2) = total_R_attL;

cc_ipsi = squeeze(mean(cc_ipsi,3));
cc_contra = squeeze(mean(cc_contra,3));

[~,first] = find(power_left.time >= 0 & power_left.time < 0.2);
[~,second] = find(power_left.time >= 0.2 & power_left.time < 0.4);
[~,third] = find(power_left.time >= 0.4 & power_left.time < 0.6);
[~,fourth] = find(power_left.time >= 0.6 & power_left.time < 0.8);
[~,fifth] = find(power_left.time >= 0.8 & power_left.time < 1);

% Auxiliary plotting variable
positions = [0.075 0.8 0.42 0.15;...
             0.55 0.8 0.42 0.15]; 

subjs = [1 2; 3 4; 5 6; 7 8; 9 10];

% Loop through pariticipants and place them on the figure. 

for j = 1:5
    
    f = figure; 
    mp = get(0, 'MonitorPositions');
    set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17.6 20]*1.75);
    
    n = 0;
    
for i = subjs(j,:)
    n = n + 1;
    ax1 = axes('Position',positions(n,:));
    ax1.PositionConstraint = 'innerposition';
    clear Y
    Y{:,1}=cc_ipsi(i,first);    
    Y{:,2}=cc_contra(i,first);
    Y{:,3}=cc_ipsi(i,second);
    Y{:,4}=cc_contra(i,second);
    Y{:,5}=cc_ipsi(i,third);
    Y{:,6}=cc_contra(i,third);
    Y{:,7}=cc_ipsi(i,fourth);
    Y{:,8}=cc_contra(i,fourth);
    Y{:,9}=cc_ipsi(i,fifth);
    Y{:,10}=cc_contra(i,fifth);
    violin(Y,'x',x1.*scaling_factor,'facecolor',...
        [color_ipsi;color_contra;color_ipsi;color_contra;color_ipsi;color_contra;...
        color_ipsi;color_contra;color_ipsi;color_contra],'edgecolor','k',...
        'facealpha',1,'bw',0.04,'mc','k','medc',[]);
    xlim([0 1].*scaling_factor)
    
    if i == 9 || i == 10
        xticks([.1 .3 .5 .7 .9].*scaling_factor); xticklabels({'0.5–0.7','0.7–0.9','0.9–1.1','1.1–1.3','1.3–1.5'})
        xlabel('Time (s)');
        axes_handle = gca;
        axes_handle.XLabel.Position = [3.75 -0.15 -1];
    else
        xticks([]); xticklabels({''})
    end
    
    if n == 1
        yticks([0.1 0.3 0.5 0.7]); yticklabels({' .1','.3','.5','.7'})
    else
        yticks([]); yticklabels({''})
    end
    ylim([0 0.7])
    if i == 10
        str = sprintf('P%s – PLV',num2str(i));
    else
        str = sprintf('P0%s – PLV',num2str(i));
    end
    ylabel(str,'fontweight','bold')
    set(ax1,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')
   
end

%% Save figure ––––––– Uncomment and edit to save to personalised location 
% str = sprintf('s0%s_s0%s.pdf',num2str(i-1),num2str(i));
% cd '/Users/martinesparzaiaizzo/Desktop/@martirene/Figures/Figs_temp/Extended data/Figure 3-1'
% exportgraphics(gcf,str,'BackgroundColor','white','ContentType','vector')

end