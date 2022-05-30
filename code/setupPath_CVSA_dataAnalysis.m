function [folders] = setupPath_CVSA_dataAnalysis()
% 
% Martín Esparza-Iaizzo (MRG group), last version 30/05/2022
%
% FieldTrip folder
    % Directorio Martin: 
    folders.Fieldtrip_folder = '/Users/martinesparzaiaizzo/Documents/MATLAB/fieldtrip-20201009';
%     folders.Fieldtrip_folder = 'D:\Martin\matlab\fieldtrip-20201009';
    addpath(folders.Fieldtrip_folder);
    ft_defaults
    fprintf('----------------------------------------\n');
    fprintf('Adding FieldTrip directory to path:\n');
    fprintf(' %s\n', folders.Fieldtrip_folder);
    
end