
%% Figures @alpha_synchronization
%
% Multisensory Research Group (MRG), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
% 
% Figure 3B. Exploratory cue-locked analysis. 
%
% Created: Thu 16 Feb 2022, 16:07
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Mon 30 May 2022, 15:32
% Last edited by: Martin Esparza-Iaizzo
% 
%
%% Load data

% load ('Exploratory Data')

%% Plot Collapsed heatmap

fontsize = 16;

[X,Y] = meshgrid((-0.5:0.002:1.5),2.^[1.65:.25:5.4]);

clear cc_ipsi cc_contra
freq_attR_R = zeros(16,1001);
freq_attR_L = zeros(16,1001);
freq_attL_R = zeros(16,1001);
freq_attL_L = zeros(16,1001);
 
for f = 1:16
    freq_attR_R(f,:) = mean(cell2mat(ss_R_attR(:,f)));
    freq_attR_L(f,:) = mean(cell2mat(ss_L_attR(:,f)));
    freq_attL_R(f,:) = mean(cell2mat(ss_R_attL(:,f)));
    freq_attL_L(f,:) = mean(cell2mat(ss_L_attL(:,f)));
end   
cc_ipsi(:,:,1) = freq_attR_R;
cc_ipsi(:,:,2) = freq_attL_L;
cc_ipsi = mean(cc_ipsi,3);
 
cc_contra(:,:,1) = freq_attR_L;
cc_contra(:,:,2) = freq_attL_R;
cc_contra = mean(cc_contra,3);
 
cc_diff = cc_contra - cc_ipsi;
cc_diff = (cc_diff - mean(mean(cc_diff)))./mean(std(cc_diff));
 
f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17.6 20]*1.75);

% ax = axes('Position',[0.1 0.625 0.8 0.075]);
% ax.PositionConstraint = 'innerposition';
% 
% ax1 = axes('Position',[0.1 0.725 0.8 0.175]);
% ax1.PositionConstraint = 'innerposition';

ax2 = axes('Position',[0.1 0.35 0.8 0.2]);
ax2.PositionConstraint = 'innerposition';
h = surf(X,Y,cc_contra - cc_ipsi);
view(2)
xlabel('Time')
ylim([3 42.2])
xlim([0 1])
ylabel('Frequency (Hz)')
xlabel('Time(s)')
% title('Difference Contra/Ipsi')
yticks([10 20 30 40]); yticklabels({'10','20','30','40'})
c = colorbar('eastoutside','Ticks',[-0.03 0.03]);
% c.Label.String = 'z score';
caxis([-0.03 0.03])
shading interp
set(h,'edgecolor','none')
xticks([0 .2 .4 .6 .8 1]); xticklabels({'0.5', '0.7','0.9','1.1','1.3','1.5'})
set(ax2,'FontSize',fontsize,'Box','off','LineWidth',1.5,'FontName','Arial','Layer','top')
grid off

%% Save figure ––––––– Uncomment and edit to save to personalised location 
% set(gcf,'Renderer','Painter')
% exportgraphics(gcf,'cue_locked.eps');

