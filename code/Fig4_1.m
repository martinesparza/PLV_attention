
%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 4A. Individual Target-locked results
%
% Created: Mon 16 Nov 2020, 10:05
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Sat 12 Feb 2022, 15:31
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

%% Load all subjects

% Add target-locked data
DataFolder = '/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/datasets/dataset2_ref2target';

% Change current folder to Subject Folder
cd(DataFolder);
fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected subject \n');

% Initialize global variable
% Rows: subjects, Col 1: PL-pre, Col 2: PL-post, Col 3: PR-pre, Col 4:
% PR-post
PLV_att_right = zeros(10,4); 
PLV_att_left = zeros(10,4);

subjects = 10;

ss_PLV_pre_ipsi = cell(1,10);
ss_PLV_pre_contra = cell(1,10);
ss_PLV_post_ipsi = cell(1,10);
ss_PLV_post_contra = cell(1,10);

%% Main loop

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
cfg.channel    = {'Fz';'O2';'FC1';'P6';'FC2';'P8';'O1';'P5';'P7';'PO4';'PO3'};
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
        tmp(1,:,:,:) = squeeze(angle(power_right.fourierspctrm(:,z,:,:))); 
        tmp(2,:,:,:) = squeeze(angle(power_right.fourierspctrm(:,j,:,:))); 
        tmp = squeeze(exp(1i*diff(tmp)));
        avg_L_attR(n,:,:,:) = tmp;
        
        % Attended Left trials 
        clear tmp
        tmp(1,:,:,:) = squeeze(angle(power_left.fourierspctrm(:,z,:,:))); 
        tmp(2,:,:,:) = squeeze(angle(power_left.fourierspctrm(:,j,:,:))); 
        tmp = squeeze(exp(1i*diff(tmp)));
        avg_L_attL(n,:,:,:) = tmp;
    end
end


% FM—PR
n = 0; % Counter for positioning
for z = idx_FM'
    for j = idx_PR'
        n = n + 1;
        % Attended Right trials
        clear tmp
        tmp(1,:,:,:) = squeeze(angle(power_right.fourierspctrm(:,z,:,:))); 
        tmp(2,:,:,:) = squeeze(angle(power_right.fourierspctrm(:,j,:,:))); 
        tmp = squeeze(exp(1i*diff(tmp)));
        avg_R_attR(n,:,:,:) = tmp;
        
        % Attended Left trials 
        clear tmp
        tmp(1,:,:,:) = squeeze(angle(power_left.fourierspctrm(:,z,:,:))); 
        tmp(2,:,:,:) = squeeze(angle(power_left.fourierspctrm(:,j,:,:))); 
        tmp = squeeze(exp(1i*diff(tmp)));
        avg_R_attL(n,:,:,:) = tmp;
    end
end

% Compute consistency over time
clear final_L_attR final_L_attL final_R_attR final_R_attL 

[~,pre] = find(power_left.time >= -0.2 & power_left.time < 0);
[~,post] = find(power_left.time >= 0.2 & power_left.time < 0.4);

final_L_attR(:,:,:,1) = abs(mean(avg_L_attR(:,:,:,pre),4));
final_L_attR(:,:,:,2) = abs(mean(avg_L_attR(:,:,:,post),4));

final_L_attL(:,:,:,1) = abs(mean(avg_L_attL(:,:,:,pre),4));
final_L_attL(:,:,:,2) = abs(mean(avg_L_attL(:,:,:,post),4));

final_R_attR(:,:,:,1) = abs(mean(avg_R_attR(:,:,:,pre),4));
final_R_attR(:,:,:,2) = abs(mean(avg_R_attR(:,:,:,post),4));

final_R_attL(:,:,:,1) = abs(mean(avg_R_attL(:,:,:,pre),4));
final_R_attL(:,:,:,2) = abs(mean(avg_R_attL(:,:,:,post),4));

% Average electrode combinations and frequencies
PLV_L_attR = squeeze(mean(mean(final_L_attR,1),3));
PLV_L_attL = squeeze(mean(mean(final_L_attL,1),3));
PLV_R_attR = squeeze(mean(mean(final_R_attR,1),3));
PLV_R_attL = squeeze(mean(mean(final_R_attL,1),3));

% Average across trials. 
PLV_att_right(i,:) = [mean(PLV_L_attR) mean(PLV_R_attR)];
PLV_att_left(i,:) = [mean(PLV_L_attL) mean(PLV_R_attL)];


ss_PLV_pre_ipsi{i} = [PLV_L_attL(:,1); PLV_R_attR(:,1)];
ss_PLV_pre_contra{i} = [PLV_L_attR(:,1); PLV_R_attL(:,1)];
ss_PLV_post_ipsi{i} = [PLV_L_attL(:,2); PLV_R_attR(:,2)];
ss_PLV_post_contra{i} = [PLV_L_attR(:,2); PLV_R_attL(:,2)];


end

%% Plotting
x1 = [-0.2,-0.1,0.1,0.2];
fontsize = 16;

color_contra = [70 70 70]./255;
color_ipsi = [170 170 170]./255;

% Auxiliary plotting variables
subjs = [1 2; 3 4; 5 6; 7 8; 9 10];
positions = [0.075 0.8 0.42 0.15;...
             0.55 0.8 0.42 0.15]; 
scaling_factor = 5.5;

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
    Y{:,1}=ss_PLV_pre_ipsi{i};    
    Y{:,2}=ss_PLV_pre_contra{i};
    Y{:,3}=ss_PLV_post_ipsi{i};
    Y{:,4}=ss_PLV_post_contra{i};

    violin(Y,'x',x1.*scaling_factor,'facecolor',...
        [color_ipsi;color_contra;color_ipsi;color_contra;color_ipsi;color_contra;...
        color_ipsi;color_contra;color_ipsi;color_contra],'edgecolor','k',...
        'facealpha',1,'bw',0.04,'mc','k','medc',[]);
    xlim([-0.3 0.3].*scaling_factor)    
    xticks([-0.15 0.15].*scaling_factor); xticklabels({'Pre-target','Post-target'})
    if n == 1
        yticks([0 0.2 0.4 0.6 0.8 1]); yticklabels({'0',' .2','.4','.6','.8','1'})
    else
        yticks([]); yticklabels({''})
    end
    ylim([0 1.15])
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
% cd '/Users/martinesparzaiaizzo/Desktop/@martirene/Figures/Figs_temp/Extended data/Figure 4-1'
% exportgraphics(gcf,str,'BackgroundColor','white','ContentType','vector')
end