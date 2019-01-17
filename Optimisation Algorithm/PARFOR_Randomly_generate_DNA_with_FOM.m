
function [Generated_DNA,Corresponding_FOM] = PARFOR_Randomly_generate_DNA_with_FOM(SelectedToolbox, SelectedFigureOfMerit, use_heralding_probability, num_modes, max_ops, output_loss, DNA_length, LowerBounds, UpperBounds, Integers, Random_initial_population_size, Truncation_initial, heralding_accepted)

% Randomly generate DNA and its corresponding FOM
% Use parallel programming

%======== Rosanna Nichols and Paul Knott, 2018 =============%
%========== https://arxiv.org/abs/1812.01032 ===============%

Generated_DNA = zeros(Random_initial_population_size,DNA_length);
Corresponding_FOM = zeros(Random_initial_population_size,1);

trunc = Truncation_initial;

parfor pop = 1:Random_initial_population_size


    % MAKE RANDOM DNA
    DNA = MakeRandomDNA_noshuffle(LowerBounds, UpperBounds, Integers);

    %% Taken directly from DNAtoFOM

    %% Useful values

    max_state_settings = max(cell2mat(SelectedToolbox{1}(:,3)));
    max_operations_settings = max(cell2mat(SelectedToolbox{2}(:,3)));
    max_measurements_settings = max(cell2mat(SelectedToolbox{3}(:,3)));

    max_2_modes = floor(num_modes/2);
    num_poss_2_mode_states = sum(cell2mat(SelectedToolbox{1}(:,2))==2);
    num_measurements = num_modes - SelectedFigureOfMerit{1}{2};
    
    %% Translate DNA to variables

    current_var = 1;

    %STATES

    % vars to deal with 2-mode states. Leave out if no 2-mode states are in the toolbox.
    if num_poss_2_mode_states > 0
        %no. of 2 mode states
        num_2_mode_states = round(DNA(current_var));
        current_var = current_var + 1;
        %2 mode state selector
        two_mode_state_selectors = round(DNA(current_var:current_var + max_2_modes - 1));
        current_var = current_var + max_2_modes;

        % assignment picker - to be scaled to pick from list of possible
        % assignments of 2-mode squeezed states across all modes
        assignment_picker = DNA(current_var);
        current_var = current_var + 1;
    else
        num_2_mode_states = 0;        
    end

    % picking single mode states for each mode (settings re-used for 2-mode states)
    single_mode_state_selector = zeros(num_modes, max_state_settings + 1);
    for i = 1:num_modes
        % each row has state selector followed by settings
        single_mode_state_selector(i,:) = DNA(current_var:current_var + max_state_settings);
        current_var = current_var + max_state_settings + 1;
    end
    single_mode_state_selector(:,1) = round(single_mode_state_selector(:,1));
    
    % OPERATIONS
    operation_selector = zeros(max_ops, max_operations_settings + 3);
    for i = 1:max_ops
        % each row has [operation selector, mode selector, 2nd mode
        % selector, settings], in this order
        operation_selector(i,:) = DNA(current_var:current_var +  max_operations_settings + 2);
        current_var = current_var +  max_operations_settings + 3;
    end
    operation_selector(:,1:3) = round(operation_selector(:,1:3));
    
    % MEASUREMENTS
    measurement_selector = zeros(num_measurements, max_measurements_settings + 1);
    for i = 1:num_measurements
        % each row has [measurement selector, settings]
        measurement_selector(i,:) = DNA(current_var:current_var +  max_measurements_settings);
        current_var = current_var +  max_measurements_settings + 1;
    end
    measurement_selector(:,1) = round(measurement_selector(:,1));
    
    %% USED TO BE Loop to produce output state and if it is not truncated enough, increase the truncation
    % NOW just keep 1 trunc, and remove tests
        
    %% STATES
    % Make initial state

    states = cell(num_modes, 1);

    % Put vacuum in all modes where two mode states will be created
    if num_2_mode_states > 0

        %find which modes will have 2-mode states
        possible_assignments = MakeAssignments(num_modes, num_2_mode_states);                           % makes an array with each row being a possible assignment for which modes get 2-mode states
        int_assignment_picker = IntegerPicker (assignment_picker, 1, size(possible_assignments, 1));    % scales the assignment picker (DNA element) by the number of possible assignments
        modes_with_two_mode_states = possible_assignments(int_assignment_picker,:);                     % picks out the select assignment from the possible assignments
        chosen_2_mode_state_assignments = reshape(modes_with_two_mode_states, 2, num_2_mode_states)';    % reshapes into an array with a row for each 2-mode state

        %assign vacuum
        for i = 1:2*num_2_mode_states
            states{ modes_with_two_mode_states(i) } = FockState(0, trunc);
        end

    else
        % if there aren't any 2-mode states in this run, the list of modes with
        % 2-mode states is empty
        modes_with_two_mode_states = [];
    end

    % put assigned single mode states in modes not assigned to 2-mode states

    % select single mode states from the selected toolbox
    single_mode_state_toolbox = SelectedToolbox{1}(cell2mat(SelectedToolbox{1}(:,2))==1,:);
    for i = find(~ismember(1:num_modes, modes_with_two_mode_states)) % loop over only modes not already assigned to 2 mode states
        % if selector = 1 or more, make the state from that row of the toolbox. If selector = 0, make the vacuum
        if single_mode_state_selector(i,1) > 0
            % pick out row of toolbox for the selected state for this mode
            this_state_toolbox = single_mode_state_toolbox(single_mode_state_selector(i,1), :);
            % find the funciton that produces the state in the 10th column
            state_producing_function = this_state_toolbox{10};
            % for each setting, find the relevant element of the DNA and scale it by its limits from the toolbox
            settings = zeros(this_state_toolbox{3},1);
            for j = 1:this_state_toolbox{3}
                if this_state_toolbox{6}(j) % tests whether this setting is an integer
                    % if an integer use IntegerPicker to scale the selector from the DNA to find an integer in the range of this setting
                    settings(j) = IntegerPicker(single_mode_state_selector(i, j+1), this_state_toolbox{5}(j,1), this_state_toolbox{5}(j,2));
                else
                    % if not an integer scale the selector to be within the range of this setting
                    settings(j) = single_mode_state_selector(i, j+1) * ( this_state_toolbox{5}(j,2) - this_state_toolbox{5}(j,1) ) + this_state_toolbox{5}(j,1);
                end
            end

            % produce the state and save it in the states cell array
            % if there are fixed settings (determined by the UI input and not altered by the algorithm), include
            % these as parameters when implementing the operation
            num_fixed_settings = this_state_toolbox{7};
            if num_fixed_settings > 0
                fixed_settings = this_state_toolbox{9};
                % Produce the state with fixed settings
                states{i} = feval(state_producing_function, settings, fixed_settings, trunc);
            else
                % Produce the state without fixed settings
                states{i} = feval(state_producing_function, settings, trunc);
            end

        else
            % make vacuum
            states{i} = FockState(0, trunc);
        end
    end

    % tensor product everything together to make the initial state (with single
    % mode states in place and |0>s where the two-mode states will be produced)
    psi = TensorProduct(states);

    % act operations to form the 2 mode states on the relevant vacuum
    % states

    % select two mode states from the selected toolbox
    two_mode_state_toolbox = SelectedToolbox{1}(cell2mat(SelectedToolbox{1}(:,2))==2,:);
    for i = 1:num_2_mode_states
        modes_to_act_on = chosen_2_mode_state_assignments(i, :);
        % pick out row of toolbox for the selected state for this pair of modes
        this_state_toolbox = two_mode_state_toolbox(two_mode_state_selectors(i), :);
        % find the funciton that produces the state in the 10th column
        state_producing_function = this_state_toolbox{10};
        % for each setting, find the relevant element of the DNA (the setting variable of the first mode in the pair) and scale it by its limits from the toolbox
        settings = zeros(this_state_toolbox{3},1);
        for j = 1:this_state_toolbox{3}
            if this_state_toolbox{6}(j) % tests whether this setting is an integer
                % if an integer use IntegerPicker to scale the selector from the DNA to find an integer in the range of this setting
                settings(j) = IntegerPicker(single_mode_state_selector(modes_to_act_on(1), j+1), this_state_toolbox{5}(j,1), this_state_toolbox{5}(j,2));
            else
                % if not an integer scale the selector to be within the range of this setting
                settings(j) = single_mode_state_selector(modes_to_act_on(1), j+1) * ( this_state_toolbox{5}(j,2) - this_state_toolbox{5}(j,1) ) + this_state_toolbox{5}(j,1);
            end
        end
        % act the operation on the initial state, targeting the relevant vacuum modes to produce the desired 2-mode state
        % if there are fixed settings (determined by the UI input and not altered by the algorithm), include
        % these as parameters when implementing the operation
        num_fixed_settings = this_state_toolbox{7};
        if num_fixed_settings > 0
            fixed_settings = this_state_toolbox{9};
            % Implement the operation with fixed settings
            psi = feval(state_producing_function, psi, settings, fixed_settings, trunc, modes_to_act_on, num_modes);
        else
            % Implement the operation without fixed settings
            psi = feval(state_producing_function, psi, settings, trunc, modes_to_act_on, num_modes);
        end
    end

    %% OPERATIONS

    % act the selected operation for each layer
    for i = 1:max_ops
        % if operation selector = 0, do nothing (identity). otherwise, act
        % that operation on the current state
        if operation_selector(i,1) > 0
            % pick out row of toolbox for the selected operation
            this_operation_toolbox = SelectedToolbox{2}(operation_selector(i,1), :);
            % number of modes the operation acts on
            num_modes_to_act_on = this_operation_toolbox{2};
            % which modes the operation acts on
            acting_on_modes = zeros(num_modes_to_act_on, 1);
            acting_on_modes(1) = operation_selector(i,2);
            if num_modes_to_act_on == 2
                % second mode to act on is added to the first mode, modulo the number of modes, so the same mode is not selected
                % twice (max. of the 2nd mode to act on variable is number of modes - 2 (and min is 0))
                acting_on_modes(2) = mod(acting_on_modes(1) + operation_selector(i,3), num_modes) + 1;
            end
            % find the funciton that does the operation in the 10th column
            operation_function = this_operation_toolbox{10};
            % for each setting, find the relevant element of the DNA and scale it by its limits from the toolbox
            settings = zeros(this_operation_toolbox{3},1);
            for j = 1:this_operation_toolbox{3}
                if this_operation_toolbox{6}(j) % tests whether this setting is an integer
                    % if an integer use IntegerPicker to scale the selector from the DNA to find an integer in the range of this setting
                    settings(j) = IntegerPicker(operation_selector(i,3+j), this_operation_toolbox{5}(j,1), this_operation_toolbox{5}(j,2));
                else
                    % if not an integer scale the selector to be within the range of this setting
                    settings(j) = operation_selector(i,3+j) * ( this_operation_toolbox{5}(j,2) - this_operation_toolbox{5}(j,1) ) + this_operation_toolbox{5}(j,1);
                end
            end
            % if there are fixed settings (determined by the UI input and not altered by the algorithm), include
            % these as parameters when implementing the operation
            num_fixed_settings = this_operation_toolbox{7};
            if num_fixed_settings > 0
                fixed_settings = this_operation_toolbox{9};
                % Implement the operation with fixed settings
                psi = feval(operation_function, psi, settings, fixed_settings, trunc, acting_on_modes, num_modes);
            else
                % Implement the operation without fixed settings
                psi = feval(operation_function, psi, settings, trunc, acting_on_modes, num_modes);
            end
        end
    end

    
    %% MEASUREMENTS

    POVMs = cell(num_measurements, 1);

    for i = 1:num_measurements
        % pick out row of toolbox for the selected measurement
        this_measurement_toolbox = SelectedToolbox{3}(measurement_selector(i,1), :);
        % find POVM producing function in 10th column
        POVM_producing_function = this_measurement_toolbox{10};
        % for each setting, find the relevant element of the DNA and scale it by its limits from the toolbox
        settings = zeros(this_measurement_toolbox{3},1);
        for j = 1:this_measurement_toolbox{3}
            if this_measurement_toolbox{6}(j) % tests whether this setting is an integer
                % if an integer use IntegerPicker to scale the selector from the DNA to find an integer in the range of this setting
                settings(j) = IntegerPicker(measurement_selector(i,1+j), this_measurement_toolbox{5}(j,1), this_measurement_toolbox{5}(j,2));
            else
                % if not an integer scale the selector to be within the range of this setting
                settings(j) = measurement_selector(i,1+j) * ( this_measurement_toolbox{5}(j,2) - this_measurement_toolbox{5}(j,1) ) + this_measurement_toolbox{5}(j,1);
            end
        end

        % if there are fixed settings (determined by the UI input and not altered by the algorithm), include these as parameters when producing the POVM
        num_fixed_settings = this_measurement_toolbox{7};
        if num_fixed_settings > 0
            fixed_settings = this_measurement_toolbox{9};
            % Produce the POVM with fixed settings
            POVMs{i} = feval(POVM_producing_function, settings, fixed_settings, trunc);
        else
            % Produce the POVM without fixed settings
            POVMs{i} = feval(POVM_producing_function, settings, trunc);
        end
    end


    
    %% PURE STATE INPUT: produce final FOM value
    
    % act measurements on state
    [output_density_operator, heralding_probability] = GeneralMeasurement_pure(psi, POVMs);
    
     % LOSS

    if output_loss > 0  % if there is any loss
    
        Kraus_cut_off = trunc;  % This could be smaller than trunc, in which case ApplyLossKraus would be faster, but less accurate
    
        output_density_operator = ApplyLossKraus(output_density_operator, output_loss, Kraus_cut_off);
        
    end
    
    % FIGURE OF MERIT
    
    % evaluate Figure of Merit function on output state
    FigureofMeritUnaugmentedValue = feval(SelectedFigureOfMerit{1}{3}, output_density_operator);    

    % augmenting function - augments according to POST-measurement nbar
    nbar_post_measurement = Findnbar(output_density_operator, 1);
    AugmentingFunctionValue = feval(SelectedFigureOfMerit{2}{5}, real(nbar_post_measurement));
    
    FigureOfMeritValue = FigureofMeritUnaugmentedValue * AugmentingFunctionValue;
    
    % Penalise using the heralding probability if this option is selected
    % NOTE: this has changed from before: we don't just multiply by the heralding probability
    if use_heralding_probability % If the option is selected
        if heralding_probability < heralding_accepted % If heralding probability too small
        
        heralding_penalty = 0;   % Currently just set FOM to zero if heralding not high enough
        FigureOfMeritValue = FigureOfMeritValue * heralding_penalty;
        
        end  
    end
    
    % put in a MINUS as ga minimises the function, not maximises
    FigureOfMeritValue = - real(FigureOfMeritValue);
    
    
    %% Save DNA and FOM
    Generated_DNA(pop,:) =  DNA';
    Corresponding_FOM(pop) = FigureOfMeritValue;
    
end

end
