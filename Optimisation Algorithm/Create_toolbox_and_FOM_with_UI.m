
%% Set up and UI

load('Full Toolbox.mat') % Load full toolbox
rng('shuffle') % Shuffle random numbers
use_integers = false; % Don't use integers

set(0,'DefaultFigureWindowStyle','normal'); % Don't dock figures
% Make UI for toolbox selection, return user selections in selected toolbox
[SelectedToolbox, num_modes, max_ops, trunc_max, output_loss]  = MakeUI (FullToolbox, ToolboxTypes, default_num_modes, default_max_ops, default_trunc, default_output_loss);
% Make UI for figure of merit selection, return user selections
[SelectedFigureOfMerit, use_heralding_probability] = MakeFigureOfMeritUI_exit (FiguresOfMerit);

%set(0,'DefaultFigureWindowStyle','docked'); % Dock all figures
% Display the settings selected
MakeSelectedUI (SelectedToolbox, ToolboxTypes, num_modes, max_ops, trunc_max, output_loss);
% Display the FoM settings selected
MakeSelectedFigureOfMeritUI (SelectedFigureOfMerit, use_heralding_probability)
