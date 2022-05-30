
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
%
%% Change path to project

ProjectFolder = 'H:\Unidades compartidas\martirene\data & code\analysis';
% Change current folder to Project Folder
cd(ProjectFolder);
fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected project \n');
% PATH SETUP: Set up data analysis pipeline (path folders, toolboxes...)
setupPath_CVSA_dataAnalysis;

%% Load all subjects 

%%% Directorio Martin:
DataFolder = 'H:\Unidades compartidas\martirene\data & code\datasets\dataset2_ref2target';

% Change current folder to Subject Folder
cd(DataFolder);
fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected subject \n');
load('SubjectInformation.mat');

% Initialize global variable
iter = 10;
subjects = 10;
% Rows: subjects, Col 1: PL-pre, Col 2: PL-post, Col 3: PR-pre, Col 4:
% PR-post
PLV_att_right = zeros(10,4); 
PLV_att_left = zeros(10,4);

PLV_att_right_tmp = zeros(10,4); 
PLV_att_left_tmp = zeros(10,4);

p_value_ipsi_contra_pre = zeros(subjects, iter);
p_value_ipsi_contra_post = zeros(subjects, iter);


for i = 1:subjects

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
%     clearvars -except data_power data_power_left data_power_right NumP SubjectInfo i...
%         DataAnalysis total_R_attR total_L_attR total_R_attL total_L_attL p_value_mat iter...
%         subjects SubjectInformation 
  


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

cfg             = [];
cfg.channel    = {'Fz';'O2';'FC1';'P6';'FC2';'P8';'O1';'P5';'P7'};
cfg.output      = 'fourier';%'powandcsd';
cfg.method      = 'wavelet';
cfg.pad         = 'maxperlen';
cfg.foi         = linspace(9.54, 14.31, 10);
cfg.width       = 5;
cfg.toi         = -1:0.002:1;
cfg.keeptrials  = 'yes';
cfg.avgoverfreq = 'yes';
power_left      = ft_freqanalysis(cfg, data_filt_left); % Attended Left trials
power_right     = ft_freqanalysis(cfg, data_filt_right); % Attended Right trials


%% Connectivity analysis

clear avg_L_attR avg_L_attL avg_R_attL avg_R_attR
% Average over frequency steps
% right = squeeze(mean(power_right.fourierspctrm,3));
% left = squeeze(mean(power_left.fourierspctrm,3));

% Combine channels
% Retrieve indices
FM = {'Fz';'FC1';'FC2'};
PL = {'P7';'P5';'PO3';'O1'};
PR = {'P8';'P6';'PO4';'O2'};

idx_FM = find(~cellfun(@isempty,regexp(power_left.label,strjoin(FM,'|'))));
idx_PL = find(~cellfun(@isempty,regexp(power_left.label,strjoin(PL,'|'))));
idx_PR = find(~cellfun(@isempty,regexp(power_left.label,strjoin(PR,'|'))));

% Calculate PLV over time. First compute phase differences between each
% FM—PL location. Then average over channels to obtain the network mean and
% then measure consistency across time in pre-target and post-target time
% windows. 

% FM—PL
n = 0; % Counter for positioning
for z = idx_FM'
    for j = idx_PL'
        n = n + 1;
        % Attended Right trials
        clear tmp
        tmp(1,:,:) = squeeze(angle(right(:,z,:))); 
        tmp(2,:,:) = squeeze(angle(right(:,j,:))); 
        tmp = squeeze(exp(1i*diff(tmp)));
        avg_L_attR(n,:,:) = tmp;
        
        % Attended Left trials 
        clear tmp
        tmp(1,:,:) = squeeze(angle(left(:,z,:))); 
        tmp(2,:,:) = squeeze(angle(left(:,j,:))); 
        tmp = squeeze(exp(1i*diff(tmp)));
        avg_L_attL(n,:,:) = tmp;
    end
end

% Here, compute consistency over time for each trial & electrode
% combination.
% And then we have a real number (trials, electrode combinations, pre/post, network vs. label, frequencies)
% Now average electrode combination and frequencies. 
% Now I have trials, network, label, and pre/post
% Now, permutations, and average ipsi/contra. 

final_L_attR = squeeze(mean(avg_L_attR,1));
final_L_attL = squeeze(mean(avg_L_attL,1));

% FM—PR
n = 0; % Counter for positioning
for z = idx_FM'
    for j = idx_PR'
        n = n + 1;
        % Attended Right trials
        clear tmp
        tmp(1,:,:) = squeeze(angle(right(:,z,:))); 
        tmp(2,:,:) = squeeze(angle(right(:,j,:))); 
        tmp = squeeze(exp(1i*diff(tmp)));
        avg_R_attR(n,:,:) = tmp;
        
        % Attended Left trials 
        clear tmp
        tmp(1,:,:) = squeeze(angle(left(:,z,:))); 
        tmp(2,:,:) = squeeze(angle(left(:,j,:))); 
        tmp = squeeze(exp(1i*diff(tmp)));
        avg_R_attL(n,:,:) = tmp;
    end
end
final_R_attR = squeeze(mean(avg_R_attR,1));
final_R_attL = squeeze(mean(avg_R_attL,1));

% Compute consistency over time
clear crossTime_att_right crossTime_att_left 

[~,pre] = find(power_left.time >= -0.2 & power_left.time < 0);
[~,post] = find(power_left.time >= 0.2 & power_left.time < 0.4);

% Attended left, FM—PL
crossTime_att_left(:,1) = abs(mean(final_L_attL(:,pre),2));
crossTime_att_left(:,2) = abs(mean(final_L_attL(:,post),2));
% Attended left, FM—PR
crossTime_att_left(:,3) = abs(mean(final_R_attL(:,pre),2));
crossTime_att_left(:,4) = abs(mean(final_R_attL(:,post),2));

% Attended right, FM—PL
crossTime_att_right(:,1) = abs(mean(final_L_attR(:,pre),2));
crossTime_att_right(:,2) = abs(mean(final_L_attR(:,post),2));
% Attended right, FM—PR
crossTime_att_right(:,3) = abs(mean(final_R_attR(:,pre),2));
crossTime_att_right(:,4) = abs(mean(final_R_attR(:,post),2));

% Average across trials
PLV_att_left(i,:) = mean(crossTime_att_left);
PLV_att_right(i,:) = mean(crossTime_att_right);

%% Null distribution
% Initialize variables outside parallel loop
fourier_left = power_left.fourierspctrm;
fourier_right = power_right.fourierspctrm;

% Begin parfor loop
tic
for it = 1:iter
    
    % Define and shuffle temporary variables
    tmp_two = [fourier_left; fourier_right];
    tmp_power_left = power_left;
    tmp_power_right = power_right;
    
    len_attL = size(fourier_left);
    len_attL = len_attL(1);
    len_attR = size(fourier_right);
    len_attR = len_attR(1);
    
    perm_idx = randperm(len_attL + len_attR);  
    tmp_two = tmp_two(perm_idx,:,:,:);
    
    tmp_attL = tmp_two(1:len_attL,:,:,:);
    tmp_attR = tmp_two(len_attL+1:end,:,:,:);
     
    tmp_power_left.fourierspctrm = tmp_attL;
    tmp_power_right.fourierspctrm = tmp_attR;
    % Connectivity analysis

    avg_L_attR = [];
    avg_L_attL = [];
    avg_R_attL = [];
    avg_R_attR = [];
    % Average over frequency steps
    right = squeeze(mean(tmp_power_right.fourierspctrm,3));
    left = squeeze(mean(tmp_power_left.fourierspctrm,3));

    % Calculate PLV over time. First compute phase differences between each
    % FM—PL location. Then average over channels to obtain the network mean and
    % then measure consistency across time in pre-target and post-target time
    % windows. 

    % FM—PL
    n = 0; % Counter for positioning
    for z = idx_FM'
        for j = idx_PL'
            n = n + 1;
            % Attended Right trials
            tmp = [];
            tmp(1,:,:) = squeeze(angle(right(:,z,:))); 
            tmp(2,:,:) = squeeze(angle(right(:,j,:))); 
            tmp = squeeze(exp(1i*diff(tmp)));
            avg_L_attR(n,:,:) = tmp;

            % Attended Left trials 
            tmp = [];
            tmp(1,:,:) = squeeze(angle(left(:,z,:))); 
            tmp(2,:,:) = squeeze(angle(left(:,j,:))); 
            tmp = squeeze(exp(1i*diff(tmp)));
            avg_L_attL(n,:,:) = tmp;
        end
    end
    final_L_attR = squeeze(mean(avg_L_attR,1));
    final_L_attL = squeeze(mean(avg_L_attL,1));

    % FM—PR
    n = 0; % Counter for positioning
    for z = idx_FM'
        for j = idx_PR'
            n = n + 1;
            % Attended Right trials
            tmp = [];
            tmp(1,:,:) = squeeze(angle(right(:,z,:))); 
            tmp(2,:,:) = squeeze(angle(right(:,j,:))); 
            tmp = squeeze(exp(1i*diff(tmp)));
            avg_R_attR(n,:,:) = tmp;

            % Attended Left trials 
            tmp = [];
            tmp(1,:,:) = squeeze(angle(left(:,z,:))); 
            tmp(2,:,:) = squeeze(angle(left(:,j,:))); 
            tmp = squeeze(exp(1i*diff(tmp)));
            avg_R_attL(n,:,:) = tmp;
        end
    end
    final_R_attR = squeeze(mean(avg_R_attR,1));
    final_R_attL = squeeze(mean(avg_R_attL,1));

    % Compute consistency over time
    crossTime_att_right = [];
    crossTime_att_left = [];

    [~,pre] = find(power_left.time >= -0.2 & power_left.time < 0);
    [~,post] = find(power_left.time >= 0.2 & power_left.time < 0.4);

    % Attended left, FM-PL
    crossTime_att_left(:,1) = abs(mean(final_L_attL(:,pre),2));
    crossTime_att_left(:,2) = abs(mean(final_L_attL(:,post),2));
    % Attended left, FM-PR
    crossTime_att_left(:,3) = abs(mean(final_R_attL(:,pre),2));
    crossTime_att_left(:,4) = abs(mean(final_R_attL(:,post),2));

    % Attended right, FM-PL
    crossTime_att_right(:,1) = abs(mean(final_L_attR(:,pre),2));
    crossTime_att_right(:,2) = abs(mean(final_L_attR(:,post),2));
    % Attended right, FM-PR
    crossTime_att_right(:,3) = abs(mean(final_R_attR(:,pre),2));
    crossTime_att_right(:,4) = abs(mean(final_R_attR(:,post),2));

    % Average across trials
    PLV_att_left_temp = mean(crossTime_att_left);
    PLV_att_right_temp = mean(crossTime_att_right);
    
    PLV_ipsi_null = mean([PLV_att_left_temp(1) PLV_att_left_temp(2);...
            PLV_att_right_temp(3) PLV_att_right_temp(4)]);
        
    PLV_contra_null = mean([PLV_att_left_temp(3) PLV_att_left_temp(4);...
            PLV_att_right_temp(1) PLV_att_right_temp(2)]);
        
    p_value_ipsi_contra_pre(i,it) = (PLV_contra_null(1) - PLV_ipsi_null(1));
    p_value_ipsi_contra_post(i,it) = (PLV_contra_null(2) - PLV_ipsi_null(2));
end
toc
end

%% Calculate p_value
p_averaged(1,:) = mean(p_value_ipsi_contra_pre,1);
p_averaged(2,:) = mean(p_value_ipsi_contra_post,1);

PLV_ipsi = mean([PLV_att_left(:,1) PLV_att_left(:,2);...
            PLV_att_right(:,3) PLV_att_right(:,4)]);
        
PLV_contra = mean([PLV_att_left(:,3) PLV_att_left(:,4);...
            PLV_att_right(:,1) PLV_att_right(:,2)]);
 
p_stat_ipsi_contra_pre = (PLV_contra(1) - PLV_ipsi(1));
p_stat_ipsi_contra_post = (PLV_contra(2) - PLV_ipsi(2));

idx_p_stat_pre = find(p_averaged(1,:) > p_stat_ipsi_contra_pre);
idx_p_stat_post = find(p_averaged(2,:) > p_stat_ipsi_contra_post);
p_values_final = [length(idx_p_stat_pre)/iter length(idx_p_stat_post)/iter];

