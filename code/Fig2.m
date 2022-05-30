%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 2. Target-locked results
%
% Created: Mon 16 Nov 2020, 10:05
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Sat 12 Feb 2022, 15:31
% Last edited by: Martin Esparza-Iaizzo

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

% Add target-locked data to path. 
DataFolder = '/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/datasets/dataset2_ref2target';

% Change current folder to Subject Folder
cd(DataFolder);
fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected subject \n');

subj = 10; 
tpoints = 751; 
  
total_R_attR = zeros(subj,tpoints);
total_L_attR = zeros(subj,tpoints);
total_R_attL = zeros(subj,tpoints);
total_L_attL = zeros(subj,tpoints);

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

%% Plotting

% Auxiliary variables. 
% Collapsed as ipsi/contra
cc_ipsi = [total_R_attR; total_L_attL];
cc_contra = [total_L_attR; total_R_attL];

% Time windows
[~,pre] = find(power_left.time >= -0.2 & power_left.time <= 0);
[~,post] = find(power_left.time >= 0.2 & power_left.time <= 0.4);

% Violin variables
violin_pre_ipsi = [mean(total_L_attL(:,pre),2) mean(total_R_attR(:,pre),2)];
violin_pre_ipsi = mean(violin_pre_ipsi,2);

violin_pre_contra = [mean(total_L_attR(:,pre),2) mean(total_R_attL(:,pre),2)];
violin_pre_contra = mean(violin_pre_contra,2);

violin_post_ipsi = [mean(total_L_attL(:,post),2) mean(total_R_attR(:,post),2)];
violin_post_ipsi = mean(violin_post_ipsi,2);

violin_post_contra = [mean(total_L_attR(:,post),2) mean(total_R_attL(:,post),2)];
violin_post_contra = mean(violin_post_contra,2);

% Spacing. 
x1 = [-0.15,-0.05,0.25,0.35];

% Coloring
color_contra = [70 70 70]./255;
color_ipsi = [170 170 170]./255;
shade_ipsi = [220 220 220]./255;

color_attR = [46, 89, 132]./255;
color_attL = [114, 165, 197]./255;
shade_attL = [202, 211, 232]./255;
fontsize = 16;


f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17.6 20]*1.75);

% A
%A1 Barplot
x_overall_attR_pre = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x(1).*scaling_factor;
x_overall_attL_pre = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x(2).*scaling_factor;
x_overall_attR_post = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x(3).*scaling_factor;
x_overall_attL_post = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x(4).*scaling_factor;


ax1 = axes('Position',[0.1 0.725 0.4 0.25]);
ax1.PositionConstraint = 'innerposition';
clear Y
Y{:,1}=mean(total_L_attL(:,pre),2);
Y{:,2}=mean(total_L_attR(:,pre),2);
Y{:,3}=mean(total_L_attL(:,post),2);
Y{:,4}=mean(total_L_attR(:,post),2);
violin(Y,'x',x.*scaling_factor,'facecolor',...
    [color_attL;color_attR;color_attL;color_attR],'edgecolor','k',...
    'facealpha',1,'bw',0.04,'mc','k','medc',[]); hold on;
for i = 1:10
    hold on;
    plot([x_overall_attR_pre(i); x_overall_attL_pre(i)],[Y{:,1}(i); Y{:,2}(i)],'-','Color', 'k','Linewidth',1.1);
    hold on;
    plot([x_overall_attR_post(i); x_overall_attL_post(i)],[Y{:,3}(i); Y{:,4}(i)],'-','Color', 'k','Linewidth',1.1);
end
hold on;
scatter(x_overall_attR_pre,Y{:,1},40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.4); hold on;
scatter(x_overall_attL_pre,Y{:,2},40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.4); hold on;
scatter(x_overall_attR_post,Y{:,3},40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.4); hold on;
scatter(x_overall_attL_post,Y{:,4},40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.4); hold on;
xlim([-0.5 0.5].*scaling_factor)
xticks([]); xticklabels({''})
yticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7]); yticklabels({'.1','.2','.3','.4','.5','.6','.7'})
ylim([0.1 0.75])
ylabel('Phase Coupling')
set(ax1,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

%A1 shade
ax2 = axes('Position',[0.1 0.575 0.4 0.12]);
ax2.PositionConstraint = 'innerposition';
stdshade(total_L_attL,0.5,'-',color_attL,shade_attL,power_left.time); hold on;
stdshade(total_L_attR,0.15,'-',color_attR,color_attR,power_left.time); 
ylim([0.25 0.45])
xlim([-0.5 0.5])
yticks([0.3 0.4]); yticklabels({'.3','.4'})
xticks([-.4 -.2 0 .2 .4]); xticklabels({'-0.4','-0.2', '0', '0.2','0.4'})
ylabel('Phase coupling')
xlabel('Time (s)')
% title('Attended vs. Location')
set(ax2,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

%A2 Barplot
ax4 = axes('Position',[0.54 0.725 0.4 0.25]);
ax4.PositionConstraint = 'innerposition';
clear Y
Y{:,1}=mean(total_R_attL(:,pre),2);
Y{:,2}=mean(total_R_attR(:,pre),2);
Y{:,3}=mean(total_R_attL(:,post),2);
Y{:,4}=mean(total_R_attR(:,post),2);
violin(Y,'x',[-0.15 -0.05 0.25 0.35].*scaling_factor,'facecolor',...
    [color_attL;color_attR;color_attL;color_attR],'edgecolor','k',...
    'facealpha',1,'bw',0.04,'mc','k','medc',[]);
for i = 1:10
    hold on;
    plot([x_overall_attR_pre(i); x_overall_attL_pre(i)],[Y{:,1}(i); Y{:,2}(i)],'-','Color', 'k','Linewidth',1.1);
    hold on;
    plot([x_overall_attR_post(i); x_overall_attL_post(i)],[Y{:,3}(i); Y{:,4}(i)],'-','Color', 'k','Linewidth',1.1);
end
hold on;
scatter(x_overall_attR_pre,Y{:,1},40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.4); hold on;
scatter(x_overall_attL_pre,Y{:,2},40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.4); hold on;
scatter(x_overall_attR_post,Y{:,3},40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.4); hold on;
scatter(x_overall_attL_post,Y{:,4},40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.4); hold on;
xlim([-0.5 0.5].*scaling_factor)
xticks([]); xticklabels({''})
yticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7]); yticklabels({'.1','.2','.3','.4','.5','.6','.7'})
ylim([0.1 0.75])
% ylabel('Phase Coupling')
set(ax4,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

%A2 shade
ax3 = axes('Position',[0.54 0.575 0.4 0.12]);
ax3.PositionConstraint = 'innerposition';
stdshade(total_R_attL,0.5,'-',color_attL,shade_attL,power_left.time); hold on;
stdshade(total_R_attR,0.15,'-',color_attR,color_attR,power_left.time); 
ylim([0.25 0.45])
xlim([-0.5 0.5])
yticks([0.3 0.4]); yticklabels({'',''})
xticks([-.4 -.2 0 .2 .4]); xticklabels({'-0.4','-0.2', '0', '0.2','0.4'})
xlabel('Time (s)')
% ylabel('Phase coupling')
% title('Attended vs. Location')
set(ax3,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')


% B2 violin
clear Y
Y{:,1}=violin_pre_ipsi;
Y{:,2}=violin_pre_contra;
Y{:,3}=violin_post_ipsi;
Y{:,4}=violin_post_contra;
scaling_factor = 5.5;

x = [-0.15*ones(10,1) -0.05*ones(10,1) 0.25*ones(10,1),0.35*ones(10,1)].*scaling_factor;
y = [violin_pre_ipsi violin_pre_contra violin_post_ipsi,violin_post_contra];

ax4 = axes('Position',[0.15 0.19 0.7 0.3]);
ax4.PositionConstraint = 'innerposition';

violin(Y,'x',[-0.15 -0.05 0.25 0.35].*scaling_factor,'facecolor',...
    [color_ipsi;color_contra;color_ipsi;color_contra],'edgecolor','k',...
    'facealpha',1,'bw',0.04,'mc','k','medc',[]); hold on; 
for i = 1:10
    hold on;
    plot([x_overall_attR_pre(i); x_overall_attL_pre(i)],[Y{:,1}(i); Y{:,2}(i)],'-','Color', 'k','Linewidth',1.5);
    hold on;
    plot([x_overall_attR_post(i); x_overall_attL_post(i)],[Y{:,3}(i); Y{:,4}(i)],'-','Color', 'k','Linewidth',1.5);
end
hold on;
scatter(x_overall_attR_pre,Y{:,1},85,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',2); hold on;
scatter(x_overall_attL_pre,Y{:,2},85,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',2); hold on;
scatter(x_overall_attR_post,Y{:,3},85,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',2); hold on;
scatter(x_overall_attL_post,Y{:,4},85,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',2); hold on;

xlim([-0.5 0.5].*scaling_factor)
xticks([]); xticklabels({''})
yticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7]); yticklabels({'.1','.2','.3','.4','.5','.6','.7'})
ylim([0.1 0.75])
ylabel('Phase Coupling')
set(ax4,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

% B
ax = axes('Position',[0.15 0.05 0.7 0.12]);
ax.PositionConstraint = 'innerposition';
stdshade(cc_contra,0.10,'-',color_contra,color_contra,power_left.time); hold on;
stdshade(cc_ipsi,0.7,'-',color_ipsi,shade_ipsi,power_left.time); 

xlabel('Time')
ylim([0.25 0.45])
xlim([-0.5 0.5])
yticks([0.3 0.4]); yticklabels({'.3','.4'})
xticks([-.4 -.2 0 .2 .4]); xticklabels({'-0.4','-0.2', '0', '0.2','0.4'})
ylabel('Phase coupling')
% title('Attended vs. Location')
set(ax,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

%% Save figure ––––––– Uncomment and edit to save to personalised location 

% cd '/Users/martinesparzaiaizzo/Desktop/@martirene/Figures/Figs_temp'
% exportgraphics(gcf,'target_locked.pdf','BackgroundColor','white','ContentType','vector')
