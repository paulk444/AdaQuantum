function [v] = beam_splitter_J(dimension)
% Matrix representation of a 50:50 beam splitter, where 'dimension' is
% the cutoff of the space. 
%
% See: Rubio Jiménez, J. (2020). Non-asymptotic quantum metrology: extracting
% maximum information from limited data. PhD thesis, University of Sussex,
% ISNI: 0000 0004 8504 6357 (arXiv:1912.02324).
%
% Jesús Rubio, PhD
% Created: 16th April 2018, University of Sussex
% Last update: 8th January 2021, University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk

%% 50:50 beam splitter
v=expm(-1i*0.5*pi*sparse(j1schwinger(dimension)));
v=sparse(v);
end