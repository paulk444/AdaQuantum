function operator = CreationOperatorGeneral(trunc, acting_on_mode, number_of_modes)

% Creates a creation operator, see https://arxiv.org/abs/1812.01032
% This single-mode operator acts on mode acting_on_mode
% The total number of modes must be specified. E.g. If we have a three mode system, and we want to act on modes 1,
% we choose: acting_on_mode = 1; number_of_modes = 3

%======== Rosanna Nichols and Paul Knott, 2018 =============%
%========== https://arxiv.org/abs/1812.01032 ===============%

    adag = CreationOperator(trunc);
    operator = 1;
    
    for i = 1:number_of_modes
        
        if i == acting_on_mode
            operator = kron(operator, adag);
        else
            operator = kron(operator, eye(trunc));
        end
        
    end

end