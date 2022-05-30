
%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 3-2. Cue-locked individual exploratory figures. 
%
% Created: Fri 22 Apr 2022, 16:37
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

%% Load data

%%% Directorio Martin:
DataFolder = '/Volumes/GoogleDrive-101271366273470520077/My Drive/@martirene';

% Change current folder to Subject Folder
cd(DataFolder);
fprintf('-----------------------------------------------------------------\n');
fprintf('Changing directories to the folder of the selected data \n');

% Load exploratory data
% load ('Exploratory_Data.mat')

%% Variables

% Load subject data
clear tmp_mat_avg tmp_mat_diff tmp_mat_percent cc_ipsi cc_contra tmp_mat_avg
[X,Y] = meshgrid((-0.5:0.002:1.5),2.^[1.65:.25:5.4]);

freq_attR_R = zeros(16,1001,10);
freq_attR_L = zeros(16,1001,10);
freq_attL_R = zeros(16,1001,10);
freq_attL_L = zeros(16,1001,10);

% Select pair of subjects
s = [1 2];

for i = [1 2]
for f = 1:16

    freq_attR_R(f,:,i) = cell2mat(ss_R_attR(s(i),f));
    freq_attL_R(f,:,i) = cell2mat(ss_L_attR(s(i),f));
    freq_attR_L(f,:,i) = cell2mat(ss_R_attL(s(i),f));
    freq_attL_L(f,:,i) = cell2mat(ss_L_attL(s(i),f));
   
end

tmp_cc_ipsi(:,:,1) = freq_attR_R(:,:,i);
tmp_cc_ipsi(:,:,2) = freq_attL_L(:,:,i);
cc_ipsi(:,:,i) = mean(tmp_cc_ipsi,3);

tmp_cc_contra(:,:,1) = freq_attR_L(:,:,i);
tmp_cc_contra(:,:,2) = freq_attL_R(:,:,i);
cc_contra(:,:,i) = mean(tmp_cc_contra,3);

tmp_mat_avg(:,:,1) = cc_ipsi(:,:,i);
tmp_mat_avg(:,:,2) = cc_contra(:,:,i);
tmp_mat_avg(:,:,i) = squeeze(mean(tmp_mat_avg,3));

tmp_mat_diff = cc_contra - cc_ipsi;

tmp_mat_percent = tmp_mat_diff./tmp_mat_avg*100; 
end

%% Plotting
fontsize = 16;

f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17.6 20]*1.75);

ax = axes('Position',[0.075 0.7 0.4 0.15]);
ax.PositionConstraint = 'innerposition';
h = surf(X,Y,tmp_mat_percent(:,:,1));
view(2)
ylim([3 42.2])
xlim([0 1])
str = sprintf('P07 – Freq. (Hz)');
ylabel(str,'fontweight','bold')
% xlabel('Time(s)')
yticks([10 20 30 40]); yticklabels({' 10','20','30','40'})
caxis([-30 30])
shading interp
set(h,'edgecolor','none')
% xticks([0 .2 .4 .6 .8 1]); xticklabels({'0.5', '0.7','0.9','1.1','1.3','1.5'})
xticks([]); xticklabels({''})
set(ax,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial','Layer','top')
grid off

ax1 = axes('Position',[0.525 0.7 0.45 0.15]);
ax1.PositionConstraint = 'innerposition';
h = surf(X,Y,tmp_mat_percent(:,:,2));
view(2)
ylim([3 42.2])
xlim([0 1])
str = sprintf('P08 – Freq. (Hz)');
ylabel(str,'fontweight','bold')
% xlabel('Time(s)')
% yticks([10 20 30 40]); yticklabels({' 10','20','30','40'})
c = colorbar('eastoutside');
c.Ticks = [-30 30] ; %Create 8 ticks from zero to 1
c.TickLabels = {'-30%','30%'};
yticks([]); yticklabels({''})
caxis([-30 30])
shading interp
set(h,'edgecolor','none')
% xticks([0 .2 .4 .6 .8 1]); xticklabels({'0.5', '0.7','0.9','1.1','1.3','1.5'})
xticks([]); xticklabels({''})
set(ax1,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial','Layer','top')
grid off

%% Save figure ––––––– Uncomment and edit to save to personalised location 
% cd('/Users/martinesparzaiaizzo/Desktop/@martirene/Figures/Figs_temp/Extended data/Figure 3-2')
% set(gcf,'Renderer','Painter')
% exportgraphics(gcf,'s07_s08.pdf','BackgroundColor','white','ContentType','vector');
