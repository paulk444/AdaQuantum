function [disp] = displacement_J(dimension,alpha)
% Matrix representation of the displacement operator, where 'alpha' is the
% the amount of displacement and 'dimension' is the cutoff of the space. 
%
% See: Rubio Jiménez, J.. Non-asymptotic quantum metrology: extracting
% maximum information from limited data. PhD thesis, University of Sussex,
% ISNI: 0000 0004 8504 6357 (arXiv:1912.02324).
%
% Jesús Rubio, PhD
% Created: 16th April 2018, University of Sussex
% Last update: 8th January 2021, University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk

%% Displacement operator
disp=expm(alpha*creation(dimension)-conj(alpha)*creation(dimension)');
disp=sparse(disp);
end

