function [creat] = creation_J(dimension)
% Matrix representation of the creation operator, where 'dimension' is
% the cutoff of the space. 
%
% See: Rubio Jiménez, J. Non-asymptotic quantum metrology: extracting
% maximum information from limited data. PhD thesis, University of Sussex,
% ISNI: 0000 0004 8504 6357 (arXiv:1912.02324).
%
% Jesús Rubio, PhD
% Created: 16th April 2018, University of Sussex
% Last update: 8th January 2021, University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk

%% Creation operator
creat=zeros(dimension,dimension);
for aa=1:dimension 
    for bb=1:dimension        
        if aa == (bb+1)        
        creat(aa,bb)=sqrt(bb);        
        else            
        creat(aa,bb)=0;       
        end
    end
end
creat=sparse(creat); 
end

