function [SelectedFigureOfMerit, use_heralding_probability] = GetInputFigureOfMerit(FiguresOfMerit, FigureOfMeritButtons, Heralding_prob_checkbox, AugmentingFunctionButtons, AugSettingEdits)
% reads in the user input from the UI and saves it in the selected toolbox
% (and other relevant variables) so these can be used when running the algorithm

    SelectedFigureOfMerit = cell(1, 2);

    % read in whether or not to use the heralding probability as a penalty
    % function using the checkbox
    use_heralding_probability = Heralding_prob_checkbox.Value;

    % Selected figure of merit
    for FoM_index = 1:length(FigureOfMeritButtons)
        % find the figure of merit radio button that has been selected
       if  FigureOfMeritButtons(FoM_index).Value
           % copy the relevant row from the figures of merit cell array
          SelectedFigureOfMerit{1} = FiguresOfMerit{1}(FoM_index,:);
       end
        
    end
    
    % Selected Augmenting Function
    for Aug_index = 1:length(AugmentingFunctionButtons)
        % find the figure of merit radio button that has been selected
        if  AugmentingFunctionButtons(Aug_index).Value
            % copy the relevant row from the figures of merit cell array
            SelectedFigureOfMerit{2} = FiguresOfMerit{2}(Aug_index,:);
            num_settings = SelectedFigureOfMerit{2}{2};
            % find the settings input in the UI
            augmenting_function_settings = zeros(num_settings,1);
            for j = 1:num_settings
                augmenting_function_settings(j) = str2double( AugSettingEdits{Aug_index}(j).String );
            end
            % save the settings to selected figure of merit (although this
            % is not actually used in DNAtoFigureOfMerit)
            SelectedFigureOfMerit{2}{4} = augmenting_function_settings;
            
            % get the augmenting-function-prodiucing-function from the
            % table and evaluate this at the selected settings to produce
            % the desired augmenting function
            AugmentingFunction = feval(SelectedFigureOfMerit{2}{5}, augmenting_function_settings);
            SelectedFigureOfMerit{2}{5} = AugmentingFunction;
            
        end
        
    end
    
    
    
end
