function [DNA_length, LowerBounds, UpperBounds, Integers] = FindDNASettings(SelectedToolbox, SelectedFigureOfMerit, num_modes, max_depth, use_integers)

% Find the DNA settings, which are needed to run main

%======== Rosanna Nichols and Paul Knott, 2018 =============%
%========== https://arxiv.org/abs/1812.01032 ===============%

    using_fixed_initial_state = strcmp(SelectedToolbox{1}(1,1), 'Two entangled pairs');

    max_state_settings = max(cell2mat(SelectedToolbox{1}(:,3)));
    max_operations_settings = max(cell2mat(SelectedToolbox{2}(:,3)));
    max_measurements_settings = max(cell2mat(SelectedToolbox{3}(:,3)));

    max_2_modes = floor(num_modes/2);
    num_poss_1_mode_states = sum(cell2mat(SelectedToolbox{1}(:,2))==1);
    num_poss_2_mode_states = sum(cell2mat(SelectedToolbox{1}(:,2))==2);
    num_poss_operations = size(SelectedToolbox{2},1);
    num_poss_measurements = size(SelectedToolbox{3},1);
    num_measurements = num_modes - SelectedFigureOfMerit{1}{2};     % Number of modes minus the number of modes the FOM acts on. All measurements are single mode

    if num_poss_2_mode_states > 0
        num_2_mode_state_vars = 2 + max_2_modes;
    else
        num_2_mode_state_vars = 0;
    end

    if using_fixed_initial_state
        DNA_states_length = 0;
    else
        DNA_states_length = num_modes * (max_state_settings + 1) + num_2_mode_state_vars;
    end
    DNA_operations_length = max_depth * (max_operations_settings + 3);
    DNA_measurements_length = num_measurements * (max_measurements_settings + 1);
    
    DNA_length =  DNA_states_length + DNA_operations_length + DNA_measurements_length;
    
    % Default settings
    LowerBounds = zeros(DNA_length,1);
    UpperBounds = ones(DNA_length,1);
    IntegerVariables = zeros(DNA_length,1);
    current_var = 1;
    
    %STATES
    
    if ~using_fixed_initial_state
        % vars to deal with 2-mode states. Leave out if no 2-mode states selected.
        if num_poss_2_mode_states > 0
            %no. of 2 mode states
            UpperBounds(current_var) = max_2_modes;
            IntegerVariables(current_var) = 1;
            current_var = current_var + 1;
            %2 mode state selector
            for i = 1:max_2_modes
                LowerBounds(current_var) = 1;
                UpperBounds(current_var) = num_poss_2_mode_states;
                IntegerVariables(current_var) = 1;
                current_var = current_var + 1;
            end
            % assignment picker - to be scaled to pick from list of possible
            % assignments of 2-mode states across all modes
            current_var = current_var + 1;
        end
        % picking single mode states for each mode (settings re-used for 2-mode
        % states)
        for i = 1:num_modes
            %state selector
            UpperBounds(current_var) = num_poss_1_mode_states;  % Lower bound = 0, where 0 is the vacuum
            IntegerVariables(current_var) = 1;
            current_var = current_var + 1;
            %settings
            for j = 1:max_state_settings
                current_var = current_var + 1;
            end
        end
    end
    
    % OPERATIONS

    % for each operation
    for i = 1:max_depth
        % operation selector
        UpperBounds(current_var) = num_poss_operations;
        IntegerVariables(current_var) = 1;
        current_var = current_var + 1;   
        % mode selector
        LowerBounds(current_var) = 1;
        UpperBounds(current_var) = num_modes;
        IntegerVariables(current_var) = 1;
        current_var = current_var + 1;   
        % second mode selector (only used for 2 mode operations)
        LowerBounds(current_var) = 0;
        UpperBounds(current_var) = num_modes - 2;
        IntegerVariables(current_var) = 1;
        current_var = current_var + 1;   
        %settings
        for j = 1:max_operations_settings
            current_var = current_var + 1;        
        end
    end

    % MEASUREMENTS

    % for each measurement
    for i = 1:num_measurements
        % measurement selector
        LowerBounds(current_var) = 1;
        UpperBounds(current_var) = num_poss_measurements;
        IntegerVariables(current_var) = 1;
        current_var = current_var + 1; 
        %settings
        for j = 1:max_measurements_settings
            current_var = current_var + 1;        
        end
    end
    
    % Check DNA is the correct length
    if DNA_length ~= current_var - 1
        disp('Error: length of DNA does not match number of variables accounted for')
    end

    % for non-integer variables move the limits very slightly away from the
    % stated limits (this was causing problems with e.g. IntegerPicker(0) = -1
    % if precisely 0 is allowed)
    for i = 1:DNA_length
        if ~IntegerVariables(i)
            LowerBounds(i) = LowerBounds(i) + 1e-6;
            UpperBounds(i) = UpperBounds(i) - 1e-6;
        end
    end
    
    
    if use_integers
        % convert list of 0s and 1s to list of variables that are integers
        Integers = find(IntegerVariables);
    else
        % if we aren't using integers, adjust the lower and upper bounds for
        % the integer variables and set the integer list to empty
        for i = 1:DNA_length
            if IntegerVariables(i)
                LowerBounds(i) = LowerBounds(i) - (0.5 - 1e-6);
                UpperBounds(i) = UpperBounds(i) + (0.5 - 1e-6);
            end
        end
        Integers = [];
    end

end