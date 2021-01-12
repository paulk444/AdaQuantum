function [squ] = squeeze_J(dimension,zeta)
% Matrix representation of the squeezing operator for a single mode,
% where 'zeta' is the squeezing parameter and 'dimension' the cutoff of
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

%% Squeezing operator
squ=expm(0.5*(conj(zeta)*(creation(dimension)')^2-zeta*(creation(dimension))^2));
squ=sparse(squ);
end

