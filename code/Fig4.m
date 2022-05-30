
%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 4A. Target-locked results
%
% Created: Mon 16 Nov 2020, 11:15
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Sat 12 Feb 2022, 21:29
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

% Data folder
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



end


%% Plotting 

x1 = [-0.15,-0.05,0.25,0.35];
x_ipsi_pre = 0.2.*(rand(1,10) - 0.5) + ones(1,10).*x1(1).*scaling_factor;
x_contra_pre = 0.2.*(rand(1,10) - 0.5) + ones(1,10).*x1(2).*scaling_factor;
x_ipsi_post = 0.2.*(rand(1,10) - 0.5) + ones(1,10).*x1(3).*scaling_factor;
x_contra_post = 0.2.*(rand(1,10) - 0.5) + ones(1,10).*x1(4).*scaling_factor;

color_contra = [70 70 70]./255;
color_ipsi = [170 170 170]./255;

fontsize = 16.5;

violin_pre_ipsi = mean([PLV_att_left(:,1) PLV_att_right(:,3)],2);
violin_pre_contra = mean([PLV_att_left(:,3) PLV_att_right(:,1)],2);
violin_post_ipsi = mean([PLV_att_left(:,2) PLV_att_right(:,4)],2);
violin_post_contra = mean([PLV_att_left(:,4) PLV_att_right(:,2)],2);

clear Y
Y{:,1}=violin_pre_ipsi; 
Y{:,2}=violin_pre_contra;
Y{:,3}=violin_post_ipsi;
Y{:,4}=violin_post_contra;
scaling_factor = 5;
x = [-0.15*ones(10,1) -0.05*ones(10,1) 0.25*ones(10,1),0.35*ones(10,1)].*scaling_factor;
y = [violin_pre_ipsi violin_pre_contra violin_post_ipsi,violin_post_contra];

f = figure; 
% mp = get(0, 'MonitorPositions');
% set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17.6 20]*1.75);
set(f,'units','centimeters','Position',[0 0 17.6 11]*2);

% Ipsi—Contra
ax1 = axes('Position',[0.075 0.6 0.275 0.35]);
ax1.PositionConstraint = 'innerposition';
violin(Y,'x',[-0.15 -0.05 0.25 0.35].*scaling_factor,'facecolor',...
    [color_ipsi;color_contra;color_ipsi;color_contra],'edgecolor','k',...
    'facealpha',1,'bw',0.04,'mc','k','medc',[]); hold on; 
for i = 1:10
    hold on;
    plot([x_ipsi_pre(i); x_contra_pre(i)],[Y{:,1}(i); Y{:,2}(i)],'-','Color', 'k','Linewidth',1.5);
    hold on;
    plot([x_ipsi_post(i); x_contra_post(i)],[Y{:,3}(i); Y{:,4}(i)],'-','Color', 'k','Linewidth',1.5);
end
hold on;
scatter(x_ipsi_pre,Y{:,1},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_contra_pre,Y{:,2},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_ipsi_post,Y{:,3},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
scatter(x_contra_post,Y{:,4},50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5); hold on;
xlim([-0.5 0.5].*scaling_factor)
xticks([-.4 -.2 0 .2 .4].*scaling_factor); xticklabels({'-0.4','-0.2', '0', '0.2','0.4'})
yticks([0.6 0.7 0.8 0.9 1]); yticklabels({'.6','.7','.8','.9','1'})
ylim([0.53 1.05])
ylabel('Phase Coupling')
xlabel('Time (s)')
set(ax1,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial','TickLength',[0.01 0.01])

% Load landscape
load('landscape.mat')
grid_accuracies = mean(global_grid_accuracies,3);
[X,Y] = meshgrid(g_values,C_values);

ax2 = axes('Position',[0.45 0.6 0.28 0.35]);
ax2.PositionConstraint = 'innerposition';
h = surf(X,Y,grid_accuracies./100);
% shading interp
set(h,'edgecolor','none')
view([120,20,20])
% colorbar
caxis([0 1])
xlim([1e-6 1e3]); %xlabel('Gamma values')
ylim([1e-6 1e3]); %ylabel('Margin values')
zlim([0 1]); zlabel('Validation accuracy')
yticks([1e-6 1e3]); yticklabels({'10^{-6}','10^3'})
xticks([1e-6]); xticklabels({'10^{-6}'})
% zticks([0 0.2 0.4 0.6 0.8 1.0]); zticklabels({'0','.2','.4','.6','.8','1'})
set(ax2,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial','XScale','log','YScale','log',...
    'YminorTick','on','XminorGrid','off','XminorTick','off')


ax3 = axes('Position',[0.79 0.74 0.18 0.15]);
ax3.PositionConstraint = 'innerposition';
H = surf(X,Y,grid_accuracies./100);
shading interp
view([120,10,20])
% colorbar
caxis([0.475 0.625])
xlim([1e-6 1e3]); %xlabel('Gamma values')
ylim([1e-6 1e3]); %ylabel('Margin values')
zlim([0.45 0.65]); %zlabel('Validation accuracy')
yticks([]); yticklabels({''})
xticks([]); xticklabels({''})
zticks([0.5 0.6]); zticklabels({'.5','.6'})
set(ax3,'FontSize',fontsize,'Box','on','LineWidth',1.5,'FontName','Arial','XScale','log','YScale','log')
% grid off

%% Save figure ––––––– Uncomment and edit to save to personalised location 

% cd '/Users/martinesparzaiaizzo/Desktop/@martirene/Figures/Figs_temp'
% set(gcf,'Renderer','Painter')
% exportgraphics(gcf,'classification_AB.eps');