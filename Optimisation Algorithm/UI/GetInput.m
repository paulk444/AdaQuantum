function [SelectedToolbox, num_modes, max_depth, trunc, output_loss] = GetInput(Checkboxes, MinMaxEdits, FixedSettingEdits, FullToolbox, NumModesEdit, MaxOpsEdit, TruncEdit, OutputLossEdit)
% reads in the user input from the UI and saves it in the selected toolbox
% (and other relevant variables) so these can be used when running the algorithm

%% MODULES & MODULE SETTINGS

    % SelectedToolbox is a cell array, very similar to FullToolbox, but
    % containing only the modules selected using the UI, with the settings
    % adjusted accordingly
    SelectedToolbox = cell(1,3);

    % type index is 1 for states, 2 for operations, 3 for measurements
    % for each of these, access the relevant UI elements and save the
    % results in the SelectedToolbox
    for type_index = 1:3

        no_modules = size(FullToolbox{type_index},1);

        for module_index = 1:no_modules

            no_settings = FullToolbox{type_index}{module_index, 3};
            no_fixed_settings = FullToolbox{type_index}{module_index, 7};
            
            %Check if checkbox for this module is ticked. If so, read in
            %settings and add to selected toolbox
            if (Checkboxes{type_index}(module_index).Value)
                
                % Read in the full information (with default values) for
                % this module from FullToolbox
                Module_data =  FullToolbox{type_index}(module_index, :);

                for settings_index = 1:no_settings
                    
                    % read in the min and max values for this setting from
                    % the edit boxes in the UI
                    min = str2double(MinMaxEdits{type_index}{module_index}(settings_index, 1).String);
                    max = str2double(MinMaxEdits{type_index}{module_index}(settings_index, 2).String);

                    % replace default settings with read in settings
                    Module_data{5}(settings_index, :) = [min, max];

                end
                
                for settings_index = 1:no_fixed_settings
                    
                    % read in value of any fixed settings (i.e. any that
                    % aren't changed by the algorithm) from the UI
                    val = str2double(FixedSettingEdits{type_index}{module_index}(settings_index).String);

                    Module_data{9}(settings_index) = val;

                end
                
                % add the information on this module as a row in the toolbox
                SelectedToolbox{type_index} = [SelectedToolbox{type_index}; Module_data];

            end

        end

    end
    
    
    %% OTHER SETTINGS
    
    num_modes = str2double(NumModesEdit.String);
    max_depth = str2double(MaxOpsEdit.String);
    trunc = str2double(TruncEdit.String);
    output_loss = str2double(OutputLossEdit.String);
    
end 