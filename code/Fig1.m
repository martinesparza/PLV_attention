
%% Figures @PLV_attention
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 1. Behavioral task and results
%
% Created: Mon 24 Jan 2022, 13:31
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Sat 12 Feb 2022, 15:31
% Last edited by: Martin Esparza-Iaizzo
% 
%
%% Behavioural data load

% Add path of data. Currently the data is private.
addpath(genpath('/Users/martinesparzaiaizzo/Desktop/@martirene/Figures/Data'))
addpath(genpath('/Volumes/GoogleDrive-101271366273470520077/Shared drives/martirene/data & code/analysis'))

% Subject mapping to avoid duplicated sessions. 
subj_idx = [1 2 3 4 5 7 8 9 10 13]; 

% Load subjects and select relevant indices
load 'Behavioral_results.mat'
data_matrix = Behavioral_table.Behavior;
data_matrix = data_matrix(subj_idx,:);

% Mean variable. Both hemifields
mean_dm = mean(data_matrix);
sem_dm = std(data_matrix)./sqrt(length(data_matrix));

%% Behavioural data plot

f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17.6 10]*1.75);

fontsize = 16;
x1 = [0.01, 0.09];
scaling_factor = 7.5;
scaling_factor1 = 7.5;


x_overall_att = 0.2.*(rand(1,10) - 0.5) + ones(1,10).*x1(1).*scaling_factor;
x_overall_unatt = 0.2.*(rand(1,10) - 0.5) + ones(1,10).*x1(2).*scaling_factor;


% Both hemifields
ax = axes('Position',[0.65 0.57 0.3 0.38]);
ax.PositionConstraint = 'innerposition';
Y{:,1}=data_matrix(:,1);
Y{:,2}=data_matrix(:,2);
violin(Y,'x',x1.*scaling_factor,'facecolor',...
    [[0.5 0.5 0.5];[0.75 0.75 0.75]],'edgecolor','k',...
    'facealpha',1,'bw',0.06,'mc','k','medc',[]);
for i = 1:10
    hold on;
    plot([x_overall_att(i); x_overall_unatt(i)],[data_matrix(i,1); data_matrix(i,2)],'-','Color', 'k','Linewidth',1.1);
end
hold on;
scatter(x_overall_att,data_matrix(:,1),50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.15); hold on;
scatter(x_overall_unatt,data_matrix(:,2),50,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.15);
ylim([0 1])
xlim([-0.1 0.2].*scaling_factor);
yticks([0 0.25 0.5 0.75 1]); yticklabels({'0','.25','.5','.75','1'})
xticks([0 0.1].*scaling_factor); xticklabels({'Attended','Unattended'})
ylabel('Proportion HITS')
set(ax,'FontSize',fontsize,'Box','on','LineWidth',1.5,'FontName','Arial')

% Right hemifield
x_overall_att = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(1).*scaling_factor;
x_overall_unatt = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(2).*scaling_factor;

ax1 = axes('Position',[0.81 0.15 0.14 0.305]);
ax1.PositionConstraint = 'innerposition';
Y{:,1}=data_matrix(:,5);
Y{:,2}=data_matrix(:,6);
violin(Y,'x',x1.*scaling_factor1,'facecolor',...
    [[0.5 0.5 0.5];[0.75 0.75 0.75]],'edgecolor','k',...
    'facealpha',1,'bw',0.07,'mc','k','medc',[]);hold on;
for i = 1:10
    hold on;
    plot([x_overall_att(i); x_overall_unatt(i)],[data_matrix(i,5); data_matrix(i,6)],'-','Color', 'k','Linewidth',1.1);
end
scatter(x_overall_att,data_matrix(:,5),40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.15); hold on;
scatter(x_overall_unatt,data_matrix(:,6),40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.15);
ylim([0 1.1])
xlim([-0.1 0.2].*scaling_factor1);
% yticks([0 0.5 1]); yticklabels({'0','.5','1'})
yticks([]); yticklabels({''})
xticks([0 0.1].*scaling_factor1); xticklabels({'Att.','Unatt.'})
set(ax1,'FontSize',fontsize,'Box','on','LineWidth',1.5,'FontName','Arial')


% Left hemifield
x_overall_att = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(1).*scaling_factor;
x_overall_unatt = 0.15.*(rand(1,10) - 0.5) + ones(1,10).*x1(2).*scaling_factor;

ax1 = axes('Position',[0.65 0.15 0.14 0.305]);
ax1.PositionConstraint = 'innerposition';
Y{:,1}=data_matrix(:,3);
Y{:,2}=data_matrix(:,4);
violin(Y,'x',x1.*scaling_factor1,'facecolor',...
    [[0.5 0.5 0.5];[0.75 0.75 0.75]],'edgecolor','k',...
    'facealpha',1,'bw',0.07,'mc','k','medc',[]);hold on;
for i = 1:10
    hold on;
    plot([x_overall_att(i); x_overall_unatt(i)],[data_matrix(i,3); data_matrix(i,4)],'-','Color', 'k','Linewidth',1.1);
end
scatter(x_overall_att,data_matrix(:,3),40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.15); hold on;
scatter(x_overall_unatt,data_matrix(:,4),40,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.15);
ylim([0 1.1])
xlim([-0.1 0.2].*scaling_factor1);
yticks([0 0.5 1]); yticklabels({'0','.5','1'})
xticks([0 0.1].*scaling_factor1); xticklabels({'Att.','Unatt.'})
ylabel('Proportion HITS')
set(ax1,'FontSize',fontsize,'Box','on','LineWidth',1.5,'FontName','Arial')

%% Save figure ––––––– Uncomment and edit to save to personalised location 

% cd '/Users/martinesparzaiaizzo/Desktop/@martirene/Figures/Figs_temp'
% exportgraphics(gcf,'behavioral.eps','BackgroundColor','white','ContentType','vector')

