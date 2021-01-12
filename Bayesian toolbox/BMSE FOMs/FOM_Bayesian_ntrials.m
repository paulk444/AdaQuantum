function [epsilon_trials] = FOM_Bayesian_ntrials(pure_state_two_mode,povm_choice,prior_width,prior_mean,mu_max,select_odd_shift,num_steps,tau_mc,dim_theta)
% MODIFIED from Jesús Rubio's code by Paul Knott. Details below. Removed 'shuffle', tic & toc, and non-error displays. Changed input from "state choice" to a two-mode
% pure state: pure_state_two_mode. Added settings select_odd_shift, dim_theta.

% Jesús Rubio Jiménez, PhD student
% University of Sussex
% 17th February 2018
% J.Rubio-Jimenez@sussex.ac.uk
% 
% mz_mse_trials(state_choice,povm_choice,prior_width,prior_mean,mu_max),
%
% where 'state_choice' is the label of the initial state of a Mach-Zehnder 
% interferometer, 'povm_choice' selects the measurement scheme, 'phase_width'
% is the width of a flat prior probability, 'phase_mean' is the middle point 
% of its domain and 'mu_max' is the maximum number of observations/trials/
% repetitions.
%
% This programme calculates the mean square error as a function of the number
% of observations/trials/repetitions.

% ADDED by Jesús Rubio on Jan 2021: further the details about the theory
% and algorithm behind this calculation can be found in
%   * Rubio, J., Knott, P., and Dunningham, J. (2018). Non-asymptotic analysis of quantum metrology protocols beyond the Cramér-Rao bound. Journal of Physics Communications, 2(1):015027
%   * Rubio Jiménez, J. (2020). Non-asymptotic quantum metrology: extracting maximum information from limited data. PhD thesis, University of Sussex, ISNI: 0000 0004 8504 6357 (arXiv:1912.02324).
% tic

% Seed for the random generator
% rng('shuffle') 

% Initial state
%initial_state=initial_probe(state_choice); % OLD
initial_state = pure_state_two_mode;

% Cutoff
op_cutoff=sqrt(length(initial_state));

% Relative phase shift
a=prior_mean-prior_width/2;
b=prior_mean+prior_width/2;
% dim_theta=1250; % Now a setting
theta=linspace(a,b,dim_theta);
%num_steps=125; % Now a setting
step=round(dim_theta/num_steps);
if step-round(step)~=0
    disp('Error: dim_theta divided by num_steps must be an integer')
    return
elseif num_steps<3
    disp('Error: the approximation of external theta integral needs three rectangles at least.')
    return
end

% Monte Carlo sample size
% tau_mc=1250;  % Now a setting

% Measurement scheme
[outcomes_space,proj_columns] = mz_povm_modified(initial_state,povm_choice,prior_width,prior_mean,select_odd_shift);

% State after the phase shift, final state and amplitudes
amplitudes=zeros(length(outcomes_space),dim_theta);
for z=1:dim_theta
    after_phase_shift=sparse(phase_shift_diff_J(op_cutoff,theta(z))*initial_state);
    for x=1:length(outcomes_space)        
        povm_element=proj_columns(:,x);
        amplitudes(x,z)=sparse(povm_element'*after_phase_shift);
    end             
end

% Likelihood function
likelihood=amplitudes.*conj(amplitudes);
% disp('The likelihood function has been created.')

% Prior probability
prior=ones(1,dim_theta);
prior=prior/trapz(theta,prior);

% Bayesian inference
epsilon_bar=0;
for index_real=1:step:dim_theta
    epsilon_n=zeros(1,mu_max); % Prellocate vector  
    epsilon_n_sum=zeros(1,mu_max);
    for times=1:tau_mc
        
        % Prior density function
        prob_temp=prior;
        for runs=1:mu_max
            
            % (Monte Carlo) Interferometric simulation
            prob_sim1=likelihood(:,index_real);
            cumulative1 = cumsum(prob_sim1); % Cumulative distribution function
            prob_rand1=rand; % Random selection
            auxiliar1=cumulative1-prob_rand1;
            
            for x=1:length(outcomes_space)
                if auxiliar1(x)>0
                    index1=x;
                    break
                end
            end
            
            % Posterior density function
            prob_temp=sparse(prob_temp.*likelihood(index1,:));
            
            temp = trapz(theta,prob_temp);
            
            if temp > 1e-16
                prob_temp = prob_temp./temp;
            else
                prob_temp=0;
            end
            
            % Mean square error
            theta_expe=trapz(theta,prob_temp.*theta);
            theta2_expe=trapz(theta,prob_temp.*theta.^2);
            epsilon_n(runs)=theta2_expe-theta_expe^2;  
        end
        
        % Monte Carlo sum
        epsilon_n_sum=epsilon_n_sum+epsilon_n;
    end
       
    % Monte Carlo approximation
    epsilon_average=epsilon_n_sum/(tau_mc);
    epsilon_bar=epsilon_bar+epsilon_average*prior(index_real)*(theta(2*step)-theta(step));
end
epsilon_trials=epsilon_bar;

%disp('The mean square error has been calculated.')
%toc

end