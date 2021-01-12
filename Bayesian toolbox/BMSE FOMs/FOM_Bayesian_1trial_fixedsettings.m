function FOM_mse_1t_final = FOM_Bayesian_1trial_fixedsettings(density_matrix_single_mode)
% This programme calculates 1 over the mean square error after 1 trial
% See [1] Rubio, J., Knott, P., and Dunningham, J. (2018). Non-asymptotic analysis of quantum metrology protocols beyond the Cramér-Rao bound. Journal of Physics Communications, 2(1):015027

% ADDED by Jesús Rubio: for the theory behind the optimal quantum bound when using
% the square error criterion with limited data (called 'bayes_bound' here), see
%   * Rubio, J., and Dunningham, J. (2019). Quantum metrology in the presence
%     of limited data. New Journal of Physics, 21 043037
%   * Rubio Jiménez, J. (2020). Non-asymptotic quantum metrology: extracting
%     maximum information from limited data. PhD thesis, University of Sussex,
%     ISNI: 0000 0004 8504 6357 (arXiv:1912.02324).

% Choose settings
prior_width = pi/12;% Width of a flat prior probability. E.g. pi/2 for squeezed states, pi for coherent states, pi/n for NOON states [1]
prior_mean = 0;     % Middle point of the piror's domain. E.g. zero

% Make the input state a pure state
pure_state_single_mode = DensityOperatorToStateVector(density_matrix_single_mode);

% Make two-mode state by single mode \otimes single mode
pure_state_two_mode = kron(pure_state_single_mode, pure_state_single_mode);

[bayes_bound] = mz_1trial_bound_evo(pure_state_two_mode,prior_width,prior_mean);

% We want to maximise FOMs...
FOM_mse_1t_final = 1/bayes_bound;

end
