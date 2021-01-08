function [outcomes,proj_columns] = mz_povm_modified(initial_state,povm_choice,phase_width,phase_mean,select_odd_shift)
% MODIFIED from Jesus Rubio's code by Paul Knott (see below). Input "state_choice" changed to a two-mode state initial_state. Removeed display options.
% Added input select_odd_shift. If this is 1, select the odd shift, otherwise select even shift.

% ADDED by Jesús Rubio on Jan 2021: see Rubio Jiménez, J. (2020). Non-asymptotic quantum metrology: extracting maximum information from limited data. PhD thesis, University of Sussex, ISNI: 0000 0004 8504 6357 (arXiv:1912.02324).

% Jesús Rubio Jiménez, PhD student
% University of Sussex
% 17th April 2018
% J.Rubio-Jimenez@sussex.ac.uk
%
% mz_povm_modified(state_choice,povm_choice,phase_width,phase_mean),
%
% where 'state_choice' labels the initial state, 'povm_choice' selects one
% of the two possible measurement strategies (see below), 'phase_width' is
% the width of the phase domain and 'phase_mean' is its centre.
%
% This programme generates the measurement outcomes (eigenvalues) and
% projective strategy (eigenvectors) of two measurement schemes:
%
%   1) Optimal measurement scheme for 1 trial using the approach developed
%   by Personick and Helstrom.
%
%   2) Known phase shift to make the strategy optimal around zero, 50:50 beam
%   splitter and photon counting.
% tic

if povm_choice==1
    %% Optimal measurement scheme for 1 trial
    [~,outcomes,proj_columns,~,~,~,~,~,~]=mz_optimal_1trial(initial_state,phase_width,phase_mean);    
elseif povm_choice==2
    %% 50:50 Beam splitter + photon counting
    
    % Cutoff
    op_cutoff=sqrt(length(initial_state));
    
    % Observable quantity (number of photons at each port)
    observable=kron(creation_J(op_cutoff)*creation_J(op_cutoff)',creation_J(op_cutoff)*creation_J(op_cutoff)');
    [proj_columns,outcomes_temp]=eig(full(observable));
    outcomes=zeros(1,length(outcomes_temp));
    for x=1:length(outcomes_temp)
        outcomes(x)=outcomes_temp(x,x);
    end
    
    % Phase shift (to make the scheme optimal around zero)
    odd_shift=kron(identity_J(op_cutoff),expm(1i*(pi/2)*creation_J(op_cutoff)*creation_J(op_cutoff)'));
    even_shift=kron(identity_J(op_cutoff),expm(1i*(pi/4)*creation_J(op_cutoff)*creation_J(op_cutoff)'));
    
    if select_odd_shift==1
        optimal_shift=odd_shift;
    else
        optimal_shift=even_shift;
    end
    
    % Effect of the 50:50 beam splitter
    proj_columns=optimal_shift'*beam_splitter_J(op_cutoff)*proj_columns;
end
%disp('The measurement scheme has been generated.')
%toc
end