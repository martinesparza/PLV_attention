%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 3. Cue-locked results. 
%
% Created: Sat 23 Apr 2022, 22:14
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Mon 30 May 2022, 15:32
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


%% Main loop
for i = 1:10
    

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
offset=2; % En segundos
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

% clearvars -except data_filt_left data_filt_right data_phase_left data_phase_right




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

% Collapsed ipsi/contra
cc_ipsi(:,:,1) = total_R_attR;
cc_ipsi(:,:,2) = total_L_attL;

cc_contra(:,:,1) = total_R_attL;
cc_contra(:,:,2) = total_R_attL;

cc_ipsi = squeeze(mean(cc_ipsi,3));
cc_contra = squeeze(mean(cc_contra,3));

% Time windows
[~,first] = find(power_left.time >= 0 & power_left.time < 0.2);
[~,second] = find(power_left.time >= 0.2 & power_left.time < 0.4);
[~,third] = find(power_left.time >= 0.4 & power_left.time < 0.6);
[~,fourth] = find(power_left.time >= 0.6 & power_left.time < 0.8);
[~,fifth] = find(power_left.time >= 0.8 & power_left.time < 1);

% Color and fontsize
color_contra = [70 70 70]./255;
color_ipsi = [170 170 170]./255;
shade_ipsi = [220 220 220]./255;

color_attL = [204, 153, 204]./255;
color_attR = [102, 153, 204]./255;
fontsize = 16;
scaling_factor = 7.5;

% Positioning
x1 = [0.07, 0.13, 0.27, 0.33, 0.47, 0.53, 0.67, 0.73, 0.87, 0.93];

% Introduce noise in the scattering. 
x_overall_ipsi_first = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(1).*scaling_factor;
x_overall_contra_first = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(2).*scaling_factor;

x_overall_ipsi_second = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(3).*scaling_factor;
x_overall_contra_second = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(4).*scaling_factor;

x_overall_ipsi_third = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(5).*scaling_factor;
x_overall_contra_third = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(6).*scaling_factor;

x_overall_ipsi_fourth = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(7).*scaling_factor;
x_overall_contra_fourth = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(8).*scaling_factor;

x_overall_ipsi_fifth = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(9).*scaling_factor;
x_overall_contra_fifth = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(10).*scaling_factor;


f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17.6 20]*1.75);

% A shade
ax = axes('Position',[0.1 0.625 0.8 0.075]);
ax.PositionConstraint = 'innerposition';
stdshade(cc_contra,0.1,'-',color_contra,color_contra,power_left.time); hold on;
stdshade(cc_ipsi,0.3,'-',color_ipsi,shade_ipsi,power_left.time); 
ylim([0.25 0.45])
xlim([0 1])
yticks([0.3 0.4]); yticklabels({'.3','.4'})
xticks([0 .2 .4 .6 .8 1]); xticklabels({'0.5', '0.7','0.9','1.1','1.3','1.5'})
% ylabel('Phase coupling')
xlabel('Time (s)')
set(ax,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

% A barplot
ax1 = axes('Position',[0.1 0.725 0.8 0.175]);
ax1.PositionConstraint = 'innerposition';
clear Y
Y{:,1}=mean(cc_ipsi(:,first),2);
Y{:,2}=mean(cc_contra(:,first),2);
Y{:,3}=mean(cc_ipsi(:,second),2);
Y{:,4}=mean(cc_contra(:,second),2);
Y{:,5}=mean(cc_ipsi(:,third),2);
Y{:,6}=mean(cc_contra(:,third),2);
Y{:,7}=mean(cc_ipsi(:,fourth),2);
Y{:,8}=mean(cc_contra(:,fourth),2);
Y{:,9}=mean(cc_ipsi(:,fifth),2);
Y{:,10}=mean(cc_contra(:,fifth),2);
violin(Y,'x',x1.*scaling_factor,'facecolor',...
    [color_ipsi;color_contra;color_ipsi;color_contra;color_ipsi;color_contra;...
    color_ipsi;color_contra;color_ipsi;color_contra],'edgecolor','k',...
    'facealpha',1,'bw',0.04,'mc','k','medc',[]); hold on,
for i = 1:10
    hold on;
    plot([x_overall_ipsi_first(i); x_overall_contra_first(i)],[Y{:,1}(i); Y{:,2}(i)],'-','Color', 'k','Linewidth',1.1);
    hold on;
    plot([x_overall_ipsi_second(i); x_overall_contra_second(i)],[Y{:,3}(i); Y{:,4}(i)],'-','Color', 'k','Linewidth',1.1);
    hold on;
    plot([x_overall_ipsi_third(i); x_overall_contra_third(i)],[Y{:,5}(i); Y{:,6}(i)],'-','Color', 'k','Linewidth',1.1);
    hold on;
    plot([x_overall_ipsi_fourth(i); x_overall_contra_fourth(i)],[Y{:,7}(i); Y{:,8}(i)],'-','Color', 'k','Linewidth',1.1);
    hold on;
    plot([x_overall_ipsi_fifth(i); x_overall_contra_fifth(i)],[Y{:,9}(i); Y{:,10}(i)],'-','Color', 'k','Linewidth',1.1);
end

hold on;
scatter(x_overall_ipsi_first,Y{:,1},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_overall_contra_first,Y{:,2},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_overall_ipsi_second,Y{:,3},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_overall_contra_second,Y{:,4},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_overall_ipsi_third,Y{:,5},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_overall_contra_third,Y{:,6},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_overall_ipsi_fourth,Y{:,7},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_overall_contra_fourth,Y{:,8},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_overall_ipsi_fifth,Y{:,9},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_overall_contra_fifth,Y{:,10},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;

xlim([0 1].*scaling_factor)
xticks([]); xticklabels({''})
yticks([0.1 0.3 0.5 0.7]); yticklabels({'.1','.3','.5','.7'})
ylim([0.1 0.75])
ylabel('Phase coupling')
set(ax1,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')


%% Save figure ––––––– Uncomment and edit to save to personalised location 
% cd '/Users/martinesparzaiaizzo/Desktop/@martirene/Figures/Figs_temp'
% exportgraphics(gcf,'cue_locked_A.pdf','BackgroundColor','white','ContentType','vector')