function [utheta] = phase_shift_diff_J(dimension, theta)
% Matrix representation of the unitary encoding of the unknown parameter
% 'theta' (difference of phase shifts), where 'dimension' is the cutoff of
% the space.
%
% See: Rubio Jiménez, J. (2020). Non-asymptotic quantum metrology: extracting
% maximum information from limited data. PhD thesis, University of Sussex,
% ISNI: 0000 0004 8504 6357 (arXiv:1912.02324).
%
% Jesús Rubio, PhD
% Created: 16th April 2018, University of Sussex
% Last update: 8th January 2021, University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk

%% Unitary encoding
utheta=expm(-1i*theta*j3schwinger(dimension));
utheta=sparse(utheta);
end

