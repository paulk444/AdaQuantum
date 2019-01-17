function [output_state, total_output_probability] = GeneralMeasurement(input_state, POVM_elements)

% Takes density matrix as input state. Produces density matrix as output state
% Takes an n-mode state and m single-mode *POVM elements* in a cell array of matrices.
% Measures modes 1 to m and gives output state of modes (m + 1) to n.

%===== Lana Mineh, Rosanna Nichols, and  Paul Knott 2018 =========%
%============= https://arxiv.org/abs/1812.01032 ==================%


trunc = length(POVM_elements{1});

num_modes = log(length(input_state))/log(trunc);
num_measured_modes = length(POVM_elements);

Id_out = speye(round(trunc^(num_modes-num_measured_modes)));

total_POVM_element = POVM_elements{1};
for i = 2:num_measured_modes
    total_POVM_element = kron(total_POVM_element,POVM_elements{i});
end

% ensure sparse
total_POVM_element=sparse(total_POVM_element);

%% Calculate measurement probability

%total_output_probability = trace(kron(total_POVM_element, Id_out) * input_state);      % Old method

Id_measured=speye(round(trunc^num_measured_modes));
tr_B_rho = sparse(round(trunc^num_measured_modes),round(trunc^num_measured_modes));
trunc_output=round(trunc^(num_modes-num_measured_modes));

for i = 0:trunc-1
   thing_tr = kron(Id_measured , FockState(i,trunc_output)');
   temp_tr=input_state * thing_tr';
   tr_B_rho = tr_B_rho + thing_tr * temp_tr;
end

total_output_probability = trace(total_POVM_element*tr_B_rho);
total_output_probability = full(total_output_probability);


%% Make output state
sum_tr2 = sparse(round(trunc^(num_modes-num_measured_modes)),round(trunc^(num_modes-num_measured_modes))); % changed from sparse(trunc,trunc)
trunc_measured_modes = round(trunc^(num_measured_modes));
for i = 0:trunc_measured_modes - 1
   
   thing = kron(FockState(i,trunc_measured_modes)' , Id_out);   
   temp1=thing * kron(total_POVM_element, Id_out);
   temp2=input_state * thing';   
   sum_tr2 = sum_tr2 + temp1 * temp2;
   
   % sum_tr2 = sum_tr2 + thing * kron(total_POVM_element, Id_out) * input_state * thing';       % Old method
   
end

if total_output_probability > 0
    output_state = sum_tr2 /total_output_probability;
else
    % if heralding prob = 0, set output state = |0>, to prevent problems with NaNs
    output_state = FockState(0, trunc) * FockState(0, trunc)';
end

end