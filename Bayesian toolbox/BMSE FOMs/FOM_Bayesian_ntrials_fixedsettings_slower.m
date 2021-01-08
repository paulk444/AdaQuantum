function FOM_mse_final = FOM_Bayesian_ntrials_fixedsettings_slower(density_matrix_single_mode)
% This programme calculates the NEGATIVE OF THE final mean square error after mu_max trial
% See [1] Rubio, J., Knott, P., and Dunningham, J. (2018). Non-asymptotic analysis of quantum metrology protocols beyond the Cram´er-Rao bound. Journal of Physics Communications, 2(1):015027

% Choose settings
povm_choice = 2;    % Selects the measurement scheme, either 1 or 2. Selecting 2 gives the BS then photon counting scheme; 1 gives the optimal 1-trail scheme
prior_width = pi/12;% Width of a flat prior probability. E.g. pi/2 for squeezed states, pi for coherent states, pi/n for NOON states [1]
prior_mean = 0;     % Middle point of the piror's domain. E.g. zero?
%
mu_max = 12;        % Maximum number of observations/trials/repetitions. E.g. in [1] the mode interesting things happen in 1:100, or 2:30 more specifically
%
select_odd_shift=0; % If this is 1, select the odd shift, otherwise select even shift

% Speed settings
tau_mc = 100;        % Monte Carlo sample size, Jesus had 1250
dim_theta = 1200;     % Dimensions of the theta space to be integrated over, Jesus had 1250
num_steps = dim_theta/10;      % Number of steps in the integration?? Jesus had 125

% Make the input state a pure state
pure_state_single_mode = DensityOperatorToStateVector(density_matrix_single_mode);

% Make two-mode state by single mode \otimes single mode
pure_state_two_mode = kron(pure_state_single_mode, pure_state_single_mode);

% Call Jesús's function (modified to take pure states as input)
[epsilon_trials] = FOM_Bayesian_ntrials(pure_state_two_mode,povm_choice,prior_width,prior_mean,mu_max,select_odd_shift,num_steps,tau_mc,dim_theta);

% Find the final mean square error
mse_final = epsilon_trials(length(epsilon_trials));

% We want to maximise FOMs. We can't take the negative as the augmenting functions make the FOM closer to zero, i.e. bigger in the case of the mse.
% So, take the log, because the closer the mse is to zero, the more negative the log(mse) becomes. log(-->0) = -infty. So the closer the mse is to
% zero, the larger the negative of the log becomes!
FOM_mse_final = 1/mse_final; % NOTE: log might not be the best method!

end
