function [output_state, total_output_probability] = GeneralMeasurement_pure_test1(pure_state, POVM_elements)
% Takes PURE STATE as input state. Produces density matrix as output state
% Takes an n-mode state and m single-mode *POVM elements* in a cell array of matrices.
% Measures modes 1 to m and gives output state of modes (m + 1) to n.

% NOTE: This has only been tested for m = n-1

% Note that a POVM E_k is contructed from measurement operators M_k as
% E_k = M_k^dag M_k. Therefore, if M_k is a projector, E_k = M_k. So
% projection operators are automatically POVMs.

trunc = length(POVM_elements{1});

num_modes = log(length(pure_state))/log(trunc);
num_measured_modes = length(POVM_elements);

Id_out = speye(round(trunc^(num_modes-num_measured_modes)));

total_POVM_element = POVM_elements{1};
for i = 2:num_measured_modes
    total_POVM_element = kron(total_POVM_element,POVM_elements{i});
end

% Ensure sparse
total_POVM_element=sparse(total_POVM_element);

% Calculate measurement probability
total_output_probability = pure_state' * kron(total_POVM_element, Id_out) * pure_state;      % Old method, adapted for pure

% Make output state
sum_tr2 = sparse(round(trunc^(num_modes-num_measured_modes)),round(trunc^(num_modes-num_measured_modes)));
trunc_measured_modes = round(trunc^(num_measured_modes));
for i = 0:trunc_measured_modes - 1
    
   temp1 = kron(FockState(i,trunc_measured_modes)' * total_POVM_element , Id_out) * pure_state;
   
   temp2 = pure_state' * kron(FockState(i,trunc_measured_modes), Id_out);
    
   sum_tr2 = sum_tr2 + temp1 * temp2;
      
end

if total_output_probability > 0
    output_state = sum_tr2 /total_output_probability;
else
    % If heralding prob = 0, set output state = |0>, to prevent problems with NaNs
    output_state = FockState(0, round(trunc^(num_modes-num_measured_modes))) * FockState(0, round(trunc^(num_modes-num_measured_modes)))';
end    

if ~PurityTest(output_state, 0.001)
    
    PurityTest(output_state, 0.001)
    
end

end
