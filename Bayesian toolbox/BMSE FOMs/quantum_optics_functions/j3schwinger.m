function [j3] = j3schwinger(dimension)
% Matrix representation of the J3 operator (Jordan-Schwinger map), where 
% 'dimension' is the cutoff of the space.
%
% See: Rubio Jiménez, J. (2020). Non-asymptotic quantum metrology: extracting
% maximum information from limited data. PhD thesis, University of Sussex,
% ISNI: 0000 0004 8504 6357 (arXiv:1912.02324).
%
% Jesús Rubio, PhD
% Created: 16th April 2018, University of Sussex
% Last update: 8th January 2021, University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk

%% J3 operator (Jordan-Schwinger map)
j3=0.5*(kron(creation_J(dimension)*creation_J(dimension)',identity_J(dimension))-kron(identity_J(dimension),creation_J(dimension)*creation_J(dimension)'));
j3=sparse(j3);
end

