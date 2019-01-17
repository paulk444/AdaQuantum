function mutationChildren = mutationPower_Paul(parents,options,GenomeLength,~,~,~,thisPopulation,mutPower)
% mutationPower_Paul. Power mutation. Adapted from MUTATIONUNIFORM (the latter is built in to Matlab)
%
%   MUTATIONCHILDREN = MUTATIONPOWER_PAUL(PARENTS,OPTIONS,GENOMELENGTH,...
%                      FITNESSFCN,STATE,THISSCORE,THISPOPULATION, ...
%                      mutPower) 
% 
% This mutation function is based on K. Deep and M. Thakur, Applied mathematics and Computation 193, 211 (2007).
%
% Creates the mutated children by each element of the DNA x updating as:
%
% If t < rand, x_new = x_old - s(x_old - x_lowerbound)
% 
% If t > rand, x_new = x_old + s(x_upperbound - x_old)
% 
% where t = (x_old - x_lowerbound) / (x_upperbound - x_lowerbound)
% 
% and s = rand^mutPower
%
% i.e. mutPower = 1 means mutation just randomly moves from x_old, mutPower = infinity means no mutation at all
% For best (expected) results try: 1 < mutPower < 30
% 
%   Example:
%     options = optimoptions('ga','MutationFcn', @mutationPower_Paul); 
%
%   This will create an options structure specifying the mutation
%   function to be used is MUTATIONPOWER_PAUL.  Since the mutPower is
%   not specified, the default value of 2 is used.
%
%     mutPower = 3;
%     options = optimoptions('ga','MutationFcn', {@mutationPower_Paul, mutPower});
%
%   This will create an options structure specifying the mutation
%   function to be used is MUTATIONPOWER_PAUL and the mutPower is
%   user specified to be 3.

% NOTE: matlabs says Genome, I say DNA... they're the same???

%% INPUTS and OUTPUTS

% parents = vector of indices that specifies which elements of the population are to be mutated ???
% options = ga options
% GenomeLength = length of the DNA ?
% FitnessFcn,state,thisScore aren't used
% thisPopulation = matrix containing all the DNA of this population. Presumably each row is a DNA ???
% mutPower = power of mutation, see above

% OUTPUT = matrix where rows (?) are DNA for the children, which have been mutated from the parents' DNAs

%% This function is used in the genetic algorithm with the line: (ga-->galincon-->stepGA, line 35 of stepGA does the mutation)

% mutateKids = feval(options.MutationFcn,  parents((1 + 2 * nXoverKids):end), options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,options.MutationFcnArgs{:});

% parents((1 + 2 * nXoverKids):end) specifies which parents are mutated (the last bunch are mutated)

% options.MutationFcnArgs{:} is where mutPower is specified. Presumably it is specified by options.MutationFcn = {@mutationPower_Paul, mutPower}
% If I want to find exactly where MutationFcnArgs is specified, try searching on command line for "MutationFcnArgs" in C:\Program Files\MATLAB\R2018a\toolbox\globaloptim\globaloptim

%================== Paul Knott, 2018 =======================%
%========== https://arxiv.org/abs/1812.01032 ===============%

%%

% Set default power if not specified
if nargin < 8 || isempty(mutPower)
    mutPower = 2; % default mutation Power = 2
end

% Range of each DNA element (upper and lower bounds)
range = options.PopInitRange;
lower = range(1,:);
upper = range(2,:);
    
% I'm assuming the population is made of doubles (it's "type"??), so I've removed this if statement: "if(strcmpi(options.PopulationType,'doubleVector'))"

mutationChildren = zeros(length(parents),GenomeLength);
for i=1:length(parents) % for each parent
    
    child = thisPopulation(parents(i),:);
    
    % Mutate each element of the DNA
    
    for j=1:GenomeLength
        
        x_old = child(j);
        x_lowerbound = lower(j); % NOTE: I am assuming that PopInitRange is a matrix, as opposed to setting the range for each DNA element the same
        x_upperbound = upper(j);
        
        t = (x_old - x_lowerbound) / (x_upperbound - x_lowerbound); % used to determine whether to mutate down or up
        
        s = rand^mutPower;  % fraction of distance to mutate by
                            % this needs to be in the loops, so that each DNA element for each child gets a different random mutation
        
        if t < rand
            
            x_new = x_old - s*(x_old - x_lowerbound); % mutate left/down
            
        else
            
            x_new = x_old + s*(x_upperbound - x_old); % mutate right/up
            
        end
        
        child(j) = x_new; % update this DNA element of child
        
    end
        
    mutationChildren(i,:) = child;
    
end







