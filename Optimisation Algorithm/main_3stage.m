
%%================= AdaQuantum Main =======================%%
%--- See User Guide for how to set up and run this code ----%

%======== Rosanna Nichols and Paul Knott, 2018 =============%
%========== https://arxiv.org/abs/1812.01032 ===============%

% Load full toolbox
load('Full Toolbox.mat')


%% Things to modify

% Load selected toolbox
load('Sample_Toolbox.mat')
% Load selected FOM
load('Sample_FOM.mat')

% Set up parallel processing
p = gcp;
ga_options.UseParallel = true;

% Options for stage 1
Truncation_initial = 30;
Random_initial_population_size = 100; % Population size of stage 1

% Options for the genetic algorithm in stage 2
Population_size_s2 = 30;   % Populuation size of the first genetic algorithm in stage 2 (s2)
Elite_s2 = 3;      % Elite count s2
Cross_s2 = 0.3;    % Crossover fraction s2
mutPower = 10;      % Mutation power for both s2 and s3 (this determines how much each element of the DNA with be mutated). The larger the power, the less the mutation.
tourn_size1 = 8;    % Tournament size s2
Generations_s2 = 10; % This should be small as this stage takes yonks due to the large pop
Trunc_s2 = 80;

% Options for the genetic algorithm in stage 3
Population_size_s3 = 20; % Populuation size s3
Elite_s3 = 3;      % Elite count s3
Cross_s3 = 0.3;    % Crossover fraction s3
tourn_size2 = 8;    % Tournament size s3
StallGen_s3 = 10;  % If there is no improvement after StallGen_s3, then the ga stops
Generations_s3 = 30;   % Total number of generations for s3

% Mutation function
ga_options.MutationFcn = {@mutationPower_Paul, mutPower};

% Alternative mutation function:
% rate = 0.2
% ga_options.MutationFcn = {@mutationPower_selection, mutPower, rate};
% mutationPower_Paul mutates every element of the DNA, whereas mutationPower_selection only mutates a fraction of
% the DNA's elements. E.g. if rate = 0.1 then only 1/10th of the DNA elements will be mutated

% Crossover function
ga_options.CrossoverFcn = @crossoversinglepoint;    % Single point crossover function
% ga_options.CrossoverFcn = @crossovertwopoint;     % Two point crossover function

% Other
heralding_accepted = 1e-3;  % If the heralding probablility for a given DNA is less than this, we penalise this DNA by set the FOM value to zero



%% Set up

tic

rng('shuffle') % Shuffle random numbers
use_integers = false; % Don't use integers

% Find the DNA settings given the toolbox selected
[DNA_length, LowerBounds, UpperBounds, Integers] = FindDNASettings(SelectedToolbox, SelectedFigureOfMerit, num_modes, max_ops, use_integers);


%% Stage 1: Randomly generate DNA

disp('Starting stage 1: randomly generate strong DNA')
datetime('now')

[Generated_DNA, Corresponding_FOM] = PARFOR_Randomly_generate_DNA_with_FOM(SelectedToolbox, SelectedFigureOfMerit, use_heralding_probability, num_modes, max_ops, output_loss, DNA_length, LowerBounds, UpperBounds, Integers, Random_initial_population_size, Truncation_initial, heralding_accepted);

% Find the best DNA for the next stage
Corresponding_FOM_ordered = sort(Corresponding_FOM);    % Put the FOMs in order
Indices_to_keep = find(Corresponding_FOM <= Corresponding_FOM_ordered(Population_size_s2),Population_size_s2);   % Get indices for the top Generated_DNA -- we keep the top (ga_options.PopulationSize), to be inputted into the first ga
DNA_for_stage2 = Generated_DNA(Indices_to_keep,:); % Select DNA to keep for next stage

Generated_DNA = 0; % Set to zero as it uses a lot of memory


%% Stage 2: Run a ga with intermediate-sized population for only a few generations

disp('Starting stage 2: ga with intermediate-sized population for only a few gens')
datetime('now')
toc

% Set options chosen above
ga_options.PopulationSize = Population_size_s2;
ga_options.EliteCount = Elite_s2;
ga_options.CrossoverFraction = Cross_s2;
ga_options.Generations = Generations_s2;

ga_options.StallGenLimit = ga_options.Generations+1; % We don't want it to exit early

ga_options.InitialPopulation = DNA_for_stage2;

ga_options.PopInitRange = [LowerBounds';UpperBounds'];
ga_options.SelectionFcn = {@selectiontournament,tourn_size1};

% Make an anonymous function to pass in the parameters that are not the DNA (ga allows only the DNA as an argument)
DNAtoFigureOfMerit_s2_Anon = @(DNA) DNAtoFigureOfMerit_newSafety(DNA, SelectedToolbox, SelectedFigureOfMerit, use_heralding_probability, num_modes, max_ops, Trunc_s2, output_loss, heralding_accepted);

% run genetic algorithm
[~, ~, ~, ~, final_population, final_scores] = ga(DNAtoFigureOfMerit_s2_Anon,DNA_length,[],[],[],[],LowerBounds,UpperBounds,[],Integers,ga_options);

% adjust figure axes to something sensible
fig = gcf;
fig.Children (5).XLim = [1 inf];
fig.Children (7).XLim = [1 inf];

% FULLY SAVE the WHOLE WORKSPACE for stage 2!
% Make new folder based on current date and time
FolderName = ['Results s2, ', datestr(now,'dd mmmm yyyy HH_MM')];
mkdir(FolderName)
% Save ENTIRE workspace!
save([FolderName, '/Workspace'])

% Select the best DNA for stage 3
Corresponding_FOM_p2_results = final_scores;
Corresponding_FOM_ordered = sort(Corresponding_FOM_p2_results);    % Put the FOMs in order
Indices_to_keep = find(Corresponding_FOM_p2_results <= Corresponding_FOM_ordered(Population_size_s3),Population_size_s3);   % Get indices for the top DNA
DNA_for_stage3 = final_population(Indices_to_keep,:); % Select DNA to keep for next stage


%% Stage 3: Run a slower but more accurate ga, with a smaller population

disp('Starting stage 3: longer search with smaller pop')
datetime('now')
toc

% Set options chosen above
ga_options.PopulationSize = Population_size_s3;
ga_options.EliteCount = Elite_s3;
ga_options.CrossoverFraction = Cross_s3;
% Stopping
ga_options.StallGenLimit = StallGen_s3;
ga_options.Generations = Generations_s3;

% More-or-less fixed options for the ga
ga_options.PopInitRange = [LowerBounds';UpperBounds'];
ga_options.SelectionFcn = {@selectiontournament, tourn_size2};

% Input initial population
ga_options.InitialPopulation = DNA_for_stage3;


% Make an anonymous function to pass in the parameters that are not the DNA (ga allows only the DNA as an argument)
DNAtoFigureOfMerit_s3_Anon = @(DNA) DNAtoFigureOfMerit_newSafety(DNA, SelectedToolbox, SelectedFigureOfMerit, use_heralding_probability, num_modes, max_ops, trunc, output_loss, heralding_accepted);

% Run genetic algorithm
[WinningDNA, ~, exitflag, ga_output, final_population, final_scores] = ga(DNAtoFigureOfMerit_s3_Anon,DNA_length,[],[],[],[],LowerBounds,UpperBounds,[],Integers,ga_options);

% Adjust figure axes to something sensible
fig = gcf;
fig.Children (5).XLim = [1 inf];
fig.Children (7).XLim = [1 inf];
title('stage 3')

% Re-run the DNA to FoM function to return the output state and heralding probability as well
[FigureOfMeritValue, output_density_operator, heralding_probability, FigureofMeritUnaugmentedValue, nbar_post_measurement] = DNAtoFigureOfMerit_s3_Anon(WinningDNA);


%% Finalse

% Make new folder based on current date and time
FolderName = ['Results s3, ', datestr(now,'dd mmmm yyyy HH_MM')];
mkdir(FolderName)
% Save ENTIRE workspace!
save([FolderName, '/Workspace'])

% Print some stuff to command window whilst recording the console in a text file
diary([FolderName, '/Useful stuff print.txt'])
DNAtoDisplay(WinningDNA, SelectedToolbox, SelectedFigureOfMerit, num_modes, max_ops);
disp('Figure of merit value: ')
disp(FigureOfMeritValue)
datetime('now')
toc
diary off


%% Close parallel pool!
%delete(gcp('nocreate'))







