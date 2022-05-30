%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Created: Mon 16 Nov 2020, 10:05
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Sat 12 Feb 2022, 15:31
% Last edited by: Martin Esparza-Iaizzo
% 
%
%% Change path to project

ProjectFolder = 'H:\Unidades compartidas\martirene\data & code\analysis';
% Change current folder to Project Folder
cd(ProjectFolder);
fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected project \n');
% PATH SETUP: Set up data analysis pipeline (path folders, toolboxes...)
setupPath_CVSA_dataAnalysis;


%% Load all subjects and their IAFs

DataFolder = 'G:\Unidades compartidas\martirene\data & code\datasets\dataset2_ref2target';

% Change current folder to Subject Folder
cd(DataFolder);
fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected subject \n');
load('SubjectInformation.mat');
SubjectInformation = table2cell(SubjectInfo);

total_R_attR = zeros(10,1001);
total_L_attR = zeros(10,1001);
total_R_attL = zeros(10,1001);
total_L_attL = zeros(10,1001);

iter = 10000;

p_value_mat = zeros(iter,6,10);
subjects = [2 4 13 14 15 18 1 3 5 12];

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
        DataAnalysis total_R_attR total_L_attR total_R_attL total_L_attL p_value_mat iter...
        subjects SubjectInformation 
  


%% PREPROCESSING OF RAW DATA & TIME-FREQUENCY ANALYSIS
% Multitaper as in MBC thesis by Irene
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

IAF = cell2mat(SubjectInformation(subjects(i),2));

cfg             = [];
cfg.channel    = {'Fz';'O2';'FC1';'P6';'FC2';'P8';'O1';'P5';'P7'};
cfg.output      = 'fourier';%'powandcsd';
cfg.method      = 'wavelet';
cfg.pad         = 'maxperlen';
cfg.foi         = linspace(IAF, IAF+5, 5);
cfg.width       = 5;
cfg.toi         = -1:0.002:1;
cfg.keeptrials  = 'yes';
cfg.avgoverfreq = 'yes';
power_left      = ft_freqanalysis(cfg, data_filt_left); % Attended Left trials
power_right     = ft_freqanalysis(cfg, data_filt_right); % Attended Right trials

% Connectivity analysis

[~,pre] = find(power_left.time >= -0.2 & power_left.time <= 0);
[~,post] = find(power_left.time >= 0.2 & power_left.time <= 0.4);

% Permutations

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

tic
p_value_PL_pre = zeros(iter,1);
p_value_PL_post = zeros(iter,1);
p_value_PR_pre = zeros(iter,1);
p_value_PR_post = zeros(iter,1);
p_value_ipsi_contra_pre = zeros(iter,1);
p_value_ipsi_contra_post = zeros(iter,1);


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

    cc_ipsi = []; cc_contra = [];
    cc_ipsi = [FM_PL_attL FM_PR_attR]; cc_ipsi = mean(cc_ipsi,2);
    cc_contra = [FM_PL_attR FM_PR_attL]; cc_contra = mean(cc_contra,2);
    
    p_value_PL_pre(j) = mean(FM_PL_attR(pre)) - mean(FM_PL_attL(pre)); % Siempre contra - ipsi
    p_value_PL_post(j) = mean(FM_PL_attR(post)) - mean(FM_PL_attL(post));

    p_value_PR_pre(j) = mean(FM_PR_attL(pre)) - mean(FM_PR_attR(pre));
    p_value_PR_post(j) = mean(FM_PR_attL(post)) - mean(FM_PR_attR(post));

    p_value_ipsi_contra_pre(j) = mean(cc_contra(pre)) - mean(cc_ipsi(pre));
    p_value_ipsi_contra_post(j) = mean(cc_contra(post)) - mean(cc_ipsi(post));
end
p_value_mat(:,1,i) = p_value_PL_pre;
p_value_mat(:,2,i) = p_value_PL_post;
p_value_mat(:,3,i) = p_value_PR_pre;
p_value_mat(:,4,i) = p_value_PR_post;
p_value_mat(:,5,i) = p_value_ipsi_contra_pre;
p_value_mat(:,6,i) = p_value_ipsi_contra_post;

toc
end






