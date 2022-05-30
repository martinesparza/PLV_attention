
%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Created: Mon 16 Nov 2020, 10:05
% Author: Martin Esparza-Iaizzo
% 
% Last edited:  Mon 30 May 2022, 15:32
% Last edited by: Martin Esparza-Iaizzo
% 
%% Change path to project

% Directorio Martin:
% ProjectFolder = '/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/analysis';
ProjectFolder = 'H:\Unidades compartidas\martirene\data & code\analysis';
% Change current folder to Project Folder
cd(ProjectFolder);
fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected project \n');
% PATH SETUP: Set up data analysis pipeline (path folders, toolboxes...)
setupPath_CVSA_dataAnalysis;


%% LOAD all subjects and their IAFs

%%% Directorio Martin:
DataFolder_1 = '/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/datasets/dataset1_ref2cue';
DataFolder_0 = '/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/datasets/dataset0_rawdata';

% Change current folder to Subject Folder
addpath(DataFolder_0);
addpath(DataFolder_1);

cd(DataFolder_0)

subj = 10; 
tpoints = 151; 
total_R_attR = zeros(subj,tpoints);
total_L_attR = zeros(subj,tpoints);
total_R_attL = zeros(subj,tpoints);
total_L_attL = zeros(subj,tpoints);

subjects = [2 4 13 14 15 18 1 3 5 12];

iter = 10000;
p_value_mat = zeros(iter,subj,10); % 10 different conditions

tic
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

fourier_left = power_left.fourierspctrm;
fourier_right = power_right.fourierspctrm;

fprintf('Begin permutations... \n');


p_value_PL_first = zeros(iter,1);
p_value_PL_second = zeros(iter,1);
p_value_PL_third = zeros(iter,1);
p_value_PL_fourth = zeros(iter,1);
p_value_PL_fifth = zeros(iter,1);

p_value_PR_first = zeros(iter,1);
p_value_PR_second = zeros(iter,1);
p_value_PR_third = zeros(iter,1);
p_value_PR_fourth = zeros(iter,1);
p_value_PR_fifth = zeros(iter,1);

tic
parfor j = 1:iter
    tmp = [fourier_left; fourier_right];
    tmp_power_left = power_left;
    tmp_power_right = power_right;
    
    len_attL = size(fourier_left);
    len_attL = len_attL(1);
    len_attR = size(fourier_right);
    len_attR = len_attR(1);
    
    perm_idx = randperm(len_attL + len_attR);  
    tmp = tmp(perm_idx,:,:,:);
    
    tmp_attL = tmp(1:len_attL,:,:,:);
    tmp_attR = tmp(len_attL+1:end,:,:,:);
     

    tmp_power_left.fourierspctrm = tmp_attL;
    tmp_power_right.fourierspctrm = tmp_attR;
    

    % Fronto-medial to Parietal left ROI
    cfg             = [];
    cfg.method      = 'plv';
    cfg.channelcmb  = {'F*','O1';'F*','P5';'F*','P7'};
    FM_PL_attL = ft_connectivityanalysis(cfg, tmp_power_left);
    FM_PL_attR = ft_connectivityanalysis(cfg, tmp_power_right);

    % Fronto-medial to Parietal right ROI
    cfg             = [];
    cfg.method      = 'plv';
    cfg.channelcmb  = {'F*','O2';'F*','P6';'F*','P8'};
    FM_PR_attL = ft_connectivityanalysis(cfg, tmp_power_left);
    FM_PR_attR = ft_connectivityanalysis(cfg, tmp_power_right);

    % Average over frequencies and channels
    FM_PL_attL = squeeze(mean(FM_PL_attL.plvspctrm,1:2));
    FM_PL_attR = squeeze(mean(FM_PL_attR.plvspctrm,1:2));
    FM_PR_attL = squeeze(mean(FM_PR_attL.plvspctrm,1:2));
    FM_PR_attR = squeeze(mean(FM_PR_attR.plvspctrm,1:2));

    
    p_value_PL_first(j) = mean(FM_PL_attR(first)) - mean(FM_PL_attL(first)); % Siempre contra - ipsi
    p_value_PL_second(j) = mean(FM_PL_attR(second)) - mean(FM_PL_attL(second));
    p_value_PL_third(j) = mean(FM_PL_attR(third)) - mean(FM_PL_attL(third));
    p_value_PL_fourth(j) = mean(FM_PL_attR(fourth)) - mean(FM_PL_attL(fourth));
    p_value_PL_fifth(j) = mean(FM_PL_attR(fifth)) - mean(FM_PL_attL(fifth));

    p_value_PR_first(j) = mean(FM_PR_attL(first)) - mean(FM_PR_attR(first));
    p_value_PR_second(j) = mean(FM_PR_attL(second)) - mean(FM_PR_attR(second));
    p_value_PR_third(j) = mean(FM_PR_attL(third)) - mean(FM_PR_attR(third));
    p_value_PR_fourth(j) = mean(FM_PR_attL(fourth)) - mean(FM_PR_attR(fourth));
    p_value_PR_fifth(j) = mean(FM_PR_attL(fifth)) - mean(FM_PR_attR(fifth));

end
p_value_mat(:,1,i) = p_value_PL_first;
p_value_mat(:,2,i) = p_value_PL_second;
p_value_mat(:,3,i) = p_value_PL_third;
p_value_mat(:,4,i) = p_value_PL_fourth;
p_value_mat(:,5,i) = p_value_PL_fifth;

p_value_mat(:,6,i) = p_value_PR_first;
p_value_mat(:,7,i) = p_value_PR_second;
p_value_mat(:,8,i) = p_value_PR_third;
p_value_mat(:,9,i) = p_value_PR_fourth;
p_value_mat(:,10,i) = p_value_PR_fifth;
toc

end



