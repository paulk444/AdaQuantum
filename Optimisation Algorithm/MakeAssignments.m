function assignments = MakeAssignments(num_modes, num_2_mode_modules)

% Forms an array where each row is a possible assignment of which modes
% have 2-mode modules. The first two numbers in the row are assigned the
% first 2-mode module, the second two numbers the second module, etc.

%======== Rosanna Nichols and Paul Knott, 2018 =============%
%========== https://arxiv.org/abs/1812.01032 ===============%

    % make a vector listing the possible modes to choose from (e.g. [1,2,3,4])
    v = zeros(1,num_modes);
    for i = 1:num_modes
        v(i) = i;        
    end
    
    % for the initial list of assignments, do v choose 2 to get all the
    % possible pairs that can be chosen. If there is only one 2-mode
    % module, this is the final assingments array
    assignments = nchoosek(v, 2);
    
    % if there is more than one 2-mode module, keep finding pairs to add to
    % the intial list of pairs in assignments
    for i = 2:num_2_mode_modules
        
        % find remainders: an array where each row is the list of modes not
        % used up by that row of the assignments array
        remainders = zeros(length(assignments),num_modes-size(assignments,2));
        for j = 1:length(assignments)
           remainders(j,:) = setdiff(v,assignments(j,:));
        end
        
        % for each row in the assignments array, find all the possible
        % pairs that can be taken from the modes not already used (remainders)
        assignments_temp = [];
        for k = 1:size(remainders,1)
            % for this row in assignments, find all the possible pairs that
            % can be chosen from the remaining modes
            next_assignments = nchoosek(remainders(k,:), 2);
            for m = 1:size(next_assignments,1)
                % add a row to the new assignments array for each of the
                % new found pairs, that combines the original row from
                % assignemnts and the new pair
                assignments_temp = [assignments_temp; assignments(k,:),next_assignments(m,:)];
            end
        end
        
        assignments = assignments_temp;
        
    end
    
end