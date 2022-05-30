
%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 5. Lat. index figures. 
%
% Created: Sat 30 Apr 2022, 12:01
% Author: Irene Vigue-Guix
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

%% LOAD DATA OF INTEREST

% Variable
idx_subj = [2 4 8 9 10 13 1 3 5 7];

%% GROUP-LEVEL
% Load mean lateralization index
load('P001_mean_dp_LI.mat');
mean_dp_LI_P01 = mean_dp_LI;
load('P002_mean_dp_LI.mat');
mean_dp_LI_P02 = mean_dp_LI;
load('P003_mean_dp_LI.mat');
mean_dp_LI_P03 = mean_dp_LI;
load('P004_mean_dp_LI.mat');
mean_dp_LI_P04 = mean_dp_LI;
load('P005_mean_dp_LI.mat');
mean_dp_LI_P05 = mean_dp_LI;
load('P006_mean_dp_LI.mat');
mean_dp_LI_P06 = mean_dp_LI;
load('P007_mean_dp_LI.mat');
mean_dp_LI_P07 = mean_dp_LI;
load('P008_mean_dp_LI.mat');
mean_dp_LI_P08 = mean_dp_LI;
load('P009_mean_dp_LI.mat');
mean_dp_LI_P09 = mean_dp_LI;
load('P0010_mean_dp_LI.mat');
mean_dp_LI_P010 = mean_dp_LI;

% Adjust 
values_AttR_rightROI = [mean_dp_LI_P01.AttR; mean_dp_LI_P02.AttR; mean_dp_LI_P03.AttR; mean_dp_LI_P04.AttR; ...
    mean_dp_LI_P05.AttR; mean_dp_LI_P06.AttR; mean_dp_LI_P07.AttR; mean_dp_LI_P08.AttR; ...
    mean_dp_LI_P09.AttR; mean_dp_LI_P010.AttR];  
values_AttL_rightROI = [mean_dp_LI_P01.AttL; mean_dp_LI_P02.AttL; mean_dp_LI_P03.AttL; mean_dp_LI_P04.AttL; ...
    mean_dp_LI_P05.AttL; mean_dp_LI_P06.AttL; mean_dp_LI_P07.AttL; mean_dp_LI_P08.AttL; ...
    mean_dp_LI_P09.AttL; mean_dp_LI_P010.AttL]; 

mean_values_AttR = mean(values_AttR_rightROI);
mean_ind_values_AttR = mean(values_AttR_rightROI,2);
mean_values_AttL = mean(values_AttL_rightROI);
mean_ind_values_AttL = mean(values_AttL_rightROI,2);
std_values_AttR = std(values_AttR_rightROI);
std_values_AttL = std(values_AttL_rightROI);
sem_values_AttR  = std(values_AttR_rightROI)/sqrt(size(values_AttR_rightROI,2));
sem_values_AttL  = std(values_AttL_rightROI)/sqrt(size(values_AttL_rightROI,2));
timevec = 0:0.02:1;

%% STATS (TTEST OVER TIME)
% Stats
C1 = [mean_dp_LI_P01.AttL; mean_dp_LI_P02.AttL; mean_dp_LI_P03.AttL; ...
    mean_dp_LI_P04.AttL; mean_dp_LI_P05.AttL; mean_dp_LI_P06.AttL; ...
    mean_dp_LI_P07.AttL; mean_dp_LI_P08.AttL; mean_dp_LI_P09.AttL; ...
    mean_dp_LI_P010.AttL]; 

C2 = [mean_dp_LI_P01.AttR; mean_dp_LI_P02.AttR; mean_dp_LI_P03.AttR; ...
    mean_dp_LI_P04.AttR; mean_dp_LI_P05.AttR; mean_dp_LI_P06.AttR; ...
    mean_dp_LI_P07.AttR; mean_dp_LI_P08.AttR; mean_dp_LI_P09.AttR; ...
    mean_dp_LI_P010.AttR]; 

% PAIRED T-TEST one-tailed (left - X less than Y)
[h,Pvalue,ci,stats] = ttest(C1,C2, 'Alpha', 0.05, 'Tail','left'); 
alpha=0.05; % set alpha level
    
%%
% Mean + SEM lateralization index
figure1 = figure('WindowState','maximized','NumberTitle','off','Color',[1 1 1]);
% right
shadedErrorBar(timevec,mean_values_AttR,sem_values_AttR,'lineProps', {'LineWidth',6,...
    'LineStyle','-', 'Color',[.5 .5 .5],'MarkerFaceColor', [.5 .5 .5]}); hold on;
plot(timevec,(mean_values_AttR+sem_values_AttR),'Color',[.5 .5 .5], 'LineWidth', 3);
plot(timevec,(mean_values_AttR-sem_values_AttR),'Color',[.5 .5 .5], 'LineWidth', 3);
% left
shadedErrorBar(timevec,mean_values_AttL,sem_values_AttL,'lineProps', {'LineWidth',6,...
    'LineStyle','-', 'Color',[.96 .76 .45],'MarkerFaceColor', [.96 .76 .45]}); hold on;
plot(timevec,(mean_values_AttL+sem_values_AttL),'Color',[.96 .76 .45], 'LineWidth', 3);
plot(timevec,(mean_values_AttL-sem_values_AttL),'Color',[.96 .76 .45], 'LineWidth', 3);
v=axis; % get current axis limits
plot(timevec,(Pvalue<=alpha)*100-100+v(1)*.95,'k.','MarkerSize',25);
% plot(timevec,(c_pvalues<=c_alpha)*100-100+v(1)*.95,'k.','MarkerSize',25);
% line([0, 0], YLim, 'LineStyle','-','LineWidth', 4, 'Color', 'k');
xlabel('Time [in s]','FontSize',20);
ylim([0 0.3]);
% title('Lateralization index (mean+-sem)');
title(sprintf('Group - Lateralization Index (mean+-sem)'));
ylabel('Lat-Index','FontSize',20);
set(gca,'FontSize',28);
    % Save
    wpfilename_fig = strcat('group_mean_sem_latindex.fig');
    wpfilename_png = strcat('group_mean_sem_latindex');
    savefig(wpfilename_fig);
    print(wpfilename_png, '-dpng');
    
  
%% STATS (MEAN LAT INDEX)
% PAIRED T-TEST one-tailed (left - X less than Y)
[h,Pvalue_1,ci,stats] = ttest(mean_ind_values_AttL,mean_ind_values_AttR,...
    'Alpha', 0.05, 'Tail','left'); 

d = Cohen_d(mean_ind_values_AttL, mean_ind_values_AttR,'paired');

%% Replication of Haegens et al, (2011)
values_AttR_rightROI_1 = [mean_dp_LI_P01.AttR_rightROI; mean_dp_LI_P02.AttR_rightROI; ...
    mean_dp_LI_P03.AttR_rightROI; mean_dp_LI_P04.AttR_rightROI;...
    mean_dp_LI_P05.AttR_rightROI; mean_dp_LI_P06.AttR_rightROI; ...
    mean_dp_LI_P07.AttR_rightROI; mean_dp_LI_P08.AttR_rightROI; ...
    mean_dp_LI_P09.AttR_rightROI; mean_dp_LI_P010.AttR_rightROI];  
values_AttL_rightROI_1 = [mean_dp_LI_P01.AttL_rightROI; mean_dp_LI_P02.AttL_rightROI; ...
    mean_dp_LI_P03.AttL_rightROI; mean_dp_LI_P04.AttL_rightROI; ...
    mean_dp_LI_P05.AttL_rightROI; mean_dp_LI_P06.AttL_rightROI; ...
    mean_dp_LI_P07.AttL_rightROI; mean_dp_LI_P08.AttL_rightROI; ...
    mean_dp_LI_P09.AttL_rightROI; mean_dp_LI_P010.AttL_rightROI]; 

values_AttR_leftROI_1 = [mean_dp_LI_P01.AttR_leftROI; mean_dp_LI_P02.AttR_leftROI; ...
    mean_dp_LI_P03.AttR_leftROI; mean_dp_LI_P04.AttR_leftROI;...
    mean_dp_LI_P05.AttR_leftROI; mean_dp_LI_P06.AttR_leftROI; ...
    mean_dp_LI_P07.AttR_leftROI; mean_dp_LI_P08.AttR_leftROI; ...
    mean_dp_LI_P09.AttR_leftROI; mean_dp_LI_P010.AttR_leftROI];  
values_AttL_leftROI_1 = [mean_dp_LI_P01.AttL_leftROI; mean_dp_LI_P02.AttL_leftROI; ...
    mean_dp_LI_P03.AttL_leftROI; mean_dp_LI_P04.AttL_leftROI; ...
    mean_dp_LI_P05.AttL_leftROI; mean_dp_LI_P06.AttL_leftROI; ...
    mean_dp_LI_P07.AttL_leftROI; mean_dp_LI_P08.AttL_leftROI; ...
    mean_dp_LI_P09.AttL_leftROI; mean_dp_LI_P010.AttL_leftROI]; 

mean_values_AttR_rightROI = mean(values_AttR_rightROI_1);
mean_ind_values_AttR_rightROI = mean(values_AttR_rightROI_1,2);
mean_values_AttL_rightROI = mean(values_AttL_rightROI_1);
mean_ind_values_AttL_rightROI = mean(values_AttL_rightROI_1,2);
std_values_AttR_rightROI = std(values_AttR_rightROI_1);
std_values_AttL_rightROI = std(values_AttL_rightROI_1);
sem_values_AttR_rightROI  = std(values_AttR_rightROI_1)/sqrt(size(values_AttR_rightROI_1,2));
sem_values_AttL_rightROI  = std(values_AttL_rightROI_1)/sqrt(size(values_AttL_rightROI_1,2));

mean_values_AttR_leftROI = mean(values_AttR_leftROI_1);
mean_ind_values_AttR_leftROI = mean(values_AttR_leftROI_1,2);
mean_values_AttL_leftROI = mean(values_AttL_leftROI_1);
mean_ind_values_AttL_leftROI = mean(values_AttL_leftROI_1,2);
std_values_AttR_leftROI = std(values_AttR_leftROI_1);
std_values_AttL_leftROI = std(values_AttL_leftROI_1);
sem_values_AttR_leftROI  = std(values_AttR_leftROI_1)/sqrt(size(values_AttR_leftROI_1,2));
sem_values_AttL_leftROI  = std(values_AttL_leftROI_1)/sqrt(size(values_AttL_leftROI_1,2));
    
%% CORRELATIONS
% Change path
% cd('J:\Unitats compartides\martirene\data & code\cosas varias\datos para irene');
cd('/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/cosas varias/datos para irene');
% Load data from @Martin
load('post_contra.mat');
load('post_ipsi.mat');
load('pre_contra.mat');
load('pre_ipsi.mat');
load('SubjectInformation.mat');

% IAF
IAF_rest_power = [SubjectInfo.IAP_Power(2) SubjectInfo.IAP_Power(4) SubjectInfo.IAP_Power(8) ...
    SubjectInfo.IAP_Power(9) SubjectInfo.IAP_Power(10) SubjectInfo.IAP_Power(13) ...
    SubjectInfo.IAP_Power(1) SubjectInfo.IAP_Power(3) SubjectInfo.IAP_Power(5) ...
    SubjectInfo.IAP_Power(7)];

% IAF_rest_peak = [SubjectInfo.IAP_Hz(2) SubjectInfo.IAP_Hz(4) SubjectInfo.IAP_Hz(8) ...
%     SubjectInfo.IAP_Hz(9) SubjectInfo.IAP_Hz(10) SubjectInfo.IAP_Hz(13) ...
%     SubjectInfo.IAP_Hz(1) SubjectInfo.IAP_Hz(3) SubjectInfo.IAP_Hz(5) ...
%     SubjectInfo.IAP_Hz(7)];

IAF_task_power = [mean(mean_dp_LI_P01.power_task); mean(mean_dp_LI_P02.power_task); ...
    mean(mean_dp_LI_P03.power_task); mean(mean_dp_LI_P04.power_task); ...
    mean(mean_dp_LI_P05.power_task); mean(mean_dp_LI_P06.power_task); ...
    mean(mean_dp_LI_P07.power_task); mean(mean_dp_LI_P08.power_task); ...
    mean(mean_dp_LI_P09.power_task); mean(mean_dp_LI_P010.power_task)];

%% Effects in connectivity
% PRE
diff_pre_ipsicontra = pre_ipsi-pre_contra;
idx_plv_effect_pre = [4 5 10]; idx_plv_noeffect_pre = [1 2 3 6 7 8 9];
% POST
diff_post_ipsicontra = post_ipsi-post_contra;
idx_plv_effect_post = [4 5 6 7 10]; idx_plv_noeffect_post = [1 2 3 8 9];
     
    
%% Correction for multiple-comparisions (FDR_BH)

% FDR_BH (funciona) --> APLICADO AL MANUSCRIPT
pvals = [pval_pre pval_AttL_pre pval_AttL_pre pval_post pval_AttL_post pval_AttL_post];
q = alpha; % 0.05
method = 'pdep'; %'dep';
report = 'yes';
[h, crit_p, adj_p]=fdr_bh(pvals,q,method,report);

%% Plotting
 
color_attR = [46, 89, 132]./255;
color_attL = [114, 165, 197]./255;
shade_attL = [202, 211, 232]./255;
fontsize = 16;
scaling_factor = 7.5;
x1 = [0, 0.1];

x_overall_attL = 0.3.*(rand(1,10) - 0.5) + ones(1,10).*x1(1).*scaling_factor;
x_overall_attR = 0.3.*(rand(1,10) - 0.5) + ones(1,10).*x1(2).*scaling_factor;

%%%%%% --------

f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17.6 20]*1.75);

ax1 = axes('Position',[0.1 0.725 0.3 0.175]);
ax1.PositionConstraint = 'innerposition';
clear Y
Y{:,1}=mean_ind_values_AttL;
Y{:,2}=mean_ind_values_AttR;
violin(Y,'x',x1.*scaling_factor,'facecolor',...
    [color_attL;color_attR],'edgecolor','k',...
    'facealpha',1,'bw',0.09,'mc','k','medc',[]); hold on;
for i = 1:10
    hold on;
    plot([x_overall_attL(i); x_overall_attR(i)],[Y{:,1}(i); Y{:,2}(i)],'-','Color', 'k','Linewidth',1.75);
end
scatter(x_overall_attL,mean_ind_values_AttL,70,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.4); hold on;
scatter(x_overall_attR,mean_ind_values_AttR,70,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.4);
ylim([-0.4 0.8]);
xlim([-0.1 0.2].*scaling_factor);
yticks([-0.4 0 0.4 0.8]); yticklabels({'-.4','0','.4','.8'})
xticks([0 0.1].*scaling_factor); xticklabels({'Att. Left','Att. Right'})
ylabel('Lat. Index')
set(ax1,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

ax2 = axes('Position',[0.475 0.725 0.475 0.175]);
ax2.PositionConstraint = 'innerposition';
v = axis;
stdshade(values_AttR_rightROI,0.25,'-',color_attR,color_attR,timevec+0.5); hold on;
stdshade(values_AttL_rightROI,0.5,'-',color_attL,shade_attL,timevec+0.5); hold on;
plot(timevec+0.5,(Pvalue<=alpha)*100-100+v(1)*.95,'k.','MarkerSize',20);
xticks([0.5 0.7 0.9 1.1 1.3 1.5]); xticklabels({'0.5', '0.7','0.9','1.1','1.3','1.5'})
yticks([0.1 0.2 0.3]); yticklabels({'.1','.2','.3'})
ylabel('Lat. Index')
xlabel('Time (s)')
ylim([0 0.3])
set(ax2,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')
%% Save figure ––––––– Uncomment and edit to save to personalised location 
% exportgraphics(gcf,'irene1.pdf','BackgroundColor','white','ContentType','vector')

%% Plotting 

marker_size = 10;

%%% -- LAT iINDEX VS. Effect PRE ipsi-contra
% AttL
ax3 = axes('Position',[0.1 0.45 0.375 0.175]);
ax3.PositionConstraint = 'innerposition';
plot(mean_ind_values_AttL(idx_plv_noeffect_pre),diff_pre_ipsicontra(idx_plv_noeffect_pre), 'o', 'LineWidth',2,'MarkerSize',marker_size,...
    'MarkerEdgeColor',color_attL,'MarkerFaceColor',color_attL); hold on;
plot(mean_ind_values_AttL(idx_plv_effect_pre),diff_pre_ipsicontra(idx_plv_effect_pre), '+', 'LineWidth',2,'MarkerSize',marker_size,...
    'MarkerEdgeColor',color_attL,'MarkerFaceColor',color_attL); hold on;
% Fitting
p1 = polyfit(mean_ind_values_AttL', diff_pre_ipsicontra,1);
f1 = polyval(p1,[-5:0.01:10]);
plot([-5:0.01:10], f1, 'LineWidth', 2, 'Color', color_attL);
% Correlation
[rho_AttL_pre,pval_AttL_pre] = corr(mean_ind_values_AttL, diff_pre_ipsicontra', 'Type', 'Pearson');

% AttR    
plot(mean_ind_values_AttR(idx_plv_noeffect_pre),diff_pre_ipsicontra(idx_plv_noeffect_pre), 'o', 'LineWidth',2,'MarkerSize',marker_size,...
    'MarkerEdgeColor',color_attR,'MarkerFaceColor',color_attR); hold on;
plot(mean_ind_values_AttR(idx_plv_effect_pre),diff_pre_ipsicontra(idx_plv_effect_pre),'+', 'LineWidth',2,'MarkerSize',marker_size,...
    'MarkerEdgeColor',color_attR,'MarkerFaceColor',color_attR); hold on;
% Fitting
p2 = polyfit(mean_ind_values_AttR', diff_pre_ipsicontra,1);
f2 = polyval(p2,[-5:0.01:10]);
plot([-5:0.01:10], f2, 'LineWidth', 2, 'Color', color_attR);
% Correlation
xlim([-0.1 0.5])
ylim([-0.13 0.15])
ylabel('Effect PRE ipsi-contra');
xlabel('Lat. Index');
yticks([-0.1 0 0.1]); yticklabels({'-.1','0','.1'});
set(ax3,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

%%% -- LAT INDEX at rest VS. Effect POST ipsi-contra
% AttL
ax4 = axes('Position',[0.575 0.45 0.375 0.175]);
ax4.PositionConstraint = 'innerposition';
plot(mean_ind_values_AttL(idx_plv_noeffect_post),diff_post_ipsicontra(idx_plv_noeffect_post), 'o', 'LineWidth',2,'MarkerSize',marker_size,...
    'MarkerEdgeColor',color_attL,'MarkerFaceColor',color_attL); hold on;
plot(mean_ind_values_AttL(idx_plv_effect_post),diff_post_ipsicontra(idx_plv_effect_post), '+', 'LineWidth',2,'MarkerSize',marker_size,...
    'MarkerEdgeColor',color_attL,'MarkerFaceColor',color_attL); hold on;
% Fitting
p1 = polyfit(mean_ind_values_AttL', diff_post_ipsicontra,1);
f1 = polyval(p1,[-5:0.01:10]);
plot([-5:0.01:10], f1, 'LineWidth', 2, 'Color', color_attL);
% Correlation
[rho_AttL_pre,pval_AttL_pre] = corr(mean_ind_values_AttL, diff_post_ipsicontra', 'Type', 'Pearson');

% AttR    
plot(mean_ind_values_AttR(idx_plv_noeffect_post),diff_post_ipsicontra(idx_plv_noeffect_post), 'o', 'LineWidth',2,'MarkerSize',marker_size,...
    'MarkerEdgeColor',color_attR,'MarkerFaceColor',color_attR); hold on;
plot(mean_ind_values_AttR(idx_plv_effect_post),diff_post_ipsicontra(idx_plv_effect_post),'+', 'LineWidth',2,'MarkerSize',marker_size,...
    'MarkerEdgeColor',color_attR,'MarkerFaceColor',color_attR); hold on;
% Fitting
p2 = polyfit(mean_ind_values_AttR', diff_post_ipsicontra,1);
f2 = polyval(p2,[-5:0.01:10]);
plot([-5:0.01:10], f2, 'LineWidth', 2, 'Color', color_attR);
% Correlation
xlim([-0.1 0.5])
ylim([-0.13 0.15])
ylabel('Effect POST ipsi-contra');
xlabel('Lat. Index');
yticks([-0.1 0 0.1]); yticklabels({'-.1','0','.1'});
set(ax4,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial')

%% Save figure ––––––– Uncomment and edit to save to personalised location 
% 
% cd '/Users/martinesparzaiaizzo/Desktop/@martirene/Figures/Figs_temp'
% exportgraphics(gcf,'irene2.pdf','BackgroundColor','white','ContentType','vector')

