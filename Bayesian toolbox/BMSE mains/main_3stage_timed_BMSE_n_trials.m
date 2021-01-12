
%%================= AdaQuantum BMSE n-trials ==============%%
%--- See User Guide for how to set up and run this code ----%

% This code is specifically designed for the FOM being the Bayesian Mean Squared Error (BMSE), for n trials (see https://arxiv.org/pdf/1812.01032.pdf)

% Stage 1: Randomly generate states, and keep states with a large QFI and with nbar = 0.5
% Stage 2: BMSE FOM, n TRIALs (nbar = 0.5). Typically this FOM is very slow, so in this stage we run an approximate but faster version of the FOM, named FOM_Bayes_faster_upto05
% Stage 3: BMSE FOM, n TRIALs (nbar = 0.5). We now use the full BMSE FOM, named FOM_Bayes_upto05. Evaluating this FOM is slow, but it should be accurate.


% Load full toolbox
load('Full Toolbox.mat')

% Load selected toolbox
load('Bayes_tool_1.mat')

%--------------------------------%
p = gcp;
ga_options.UseParallel = true;
%--------------------------------%

%% Safety stuff
heralding_accepted = 1e-3;

%% Timing
tic
Time_p1_hrs = 0.01
Time_p2_hrs = 0.1
Total_time_hrs = 0.2

%% Options for stage 1
Truncation_initial = 80

%% Options for both ga's
ga_options.MutationFcn = {@mutationPower_Paul, 10}; % default 10 ?
% ga_options.MutationFcn = {@mutationPower_selection, 10, 0.05}; % default 10, 0.05 ?

% ga_options.CrossoverFcn = @crossoversinglepoint;
% ga_options.CrossoverFcn = @crossovertwopoint;

%% Options for stage 2: intermediate pop, but still do reasonably local search
Population_size_s2 = 6;
Cross_s2 = 0.3;
tourn_size1 = 8;
pop_chunk_size = Population_size_s2;
%mutPower1 = 40; % Chosen above
Elite_s2 = 1;
Generations_s2 = 20; % Set to 20, as it's running for a fixed time anyway
Trunc_s2 = 80;


%% Options for stage 3: focussed more on local search
Population_size_s3 = 2;
%
Elite_s3 = 1;
Cross_s3 = Cross_s2;
%mutPower2 = mutPower1; % Chosen above
tourn_size2 = tourn_size1;
% Other options to modify
StallGen_s3 = 50;
Generations_s3 = 300;



%% Set up and UI

rng('shuffle') % Shuffle random numbers
use_integers = false; % Don't use integers


% Find the DNA settings given the toolbox selected
[DNA_length, LowerBounds, UpperBounds, Integers] = FindDNASettings(SelectedToolbox, SelectedFigureOfMerit, num_modes, max_ops, use_integers);


%% Stage 1: Randomly generate DNA

% Load selected FOM for QFI
load('FOM_QFI_upto05.mat')

disp('Starting stage 1: randomly generate strong DNA')
datetime('now')

[Generated_DNA, Corresponding_FOM] = Randomly_generate_DNA_with_FOM_timed(SelectedToolbox, SelectedFigureOfMerit, use_heralding_probability, num_modes, max_ops, output_loss, DNA_length, LowerBounds, UpperBounds, Integers, Truncation_initial, heralding_accepted, Population_size_s2, Time_p1_hrs, pop_chunk_size);

% Find the best DNA for the next stage
Corresponding_FOM_ordered = sort(Corresponding_FOM);    % Put the FOMs in order
Indices_to_keep = find(Corresponding_FOM <= Corresponding_FOM_ordered(Population_size_s2),Population_size_s2);   % Get indices for the top Generated_DNA -- we keep the top (ga_options.PopulationSize), to be inputted into the first ga
DNA_for_stage2 = Generated_DNA(Indices_to_keep,:); % Select DNA to keep for next stage

% Save length of DNA generated, then set to zero
Random_pop_generated = size(Generated_DNA,1);
Generated_DNA = 0;

%% Stage 2: Run a ga with intermediate-sized population for only a few gens

% Load selected FOM for Bayes FASTER
load('FOM_Bayes_faster_upto05.mat')

disp('Starting stage 2: ga with intermediate-sized population for only a few gens')
datetime('now')
toc

% Set options chosen above
ga_options.PopulationSize = Population_size_s2;
ga_options.EliteCount = Elite_s2;
ga_options.CrossoverFraction = Cross_s2;
ga_options.Generations = Generations_s2;

ga_options.StallGenLimit = ga_options.Generations+1; % We don't want it to exit early
% Don't exceed time limit:
ga_options.TimeLimit = Time_p2_hrs*3600;

ga_options.InitialPopulation = DNA_for_stage2;

ga_options.PopInitRange = [LowerBounds';UpperBounds'];
ga_options.SelectionFcn = {@selectiontournament,tourn_size1};


% Make an anonymous function to pass in the parameters that are not the DNA (ga allows only the DNA as an argument)
DNAtoFigureOfMerit_s2_Anon = @(DNA) DNAtoFigureOfMerit_newSafety(DNA, SelectedToolbox, SelectedFigureOfMerit, use_heralding_probability, num_modes, max_ops, Trunc_s2, output_loss, heralding_accepted);

% run genetic algorithm
[WinningDNA, ~, exitflag, ga_output, final_population, final_scores] = ga(DNAtoFigureOfMerit_s2_Anon,DNA_length,[],[],[],[],LowerBounds,UpperBounds,[],Integers,ga_options);


% % adjust figure axes to something sensible
% fig = gcf;
% fig.Children (5).XLim = [1 inf];
% fig.Children (7).XLim = [1 inf];


% Save fig, to be reloaded later
title('Stage 2')
Fig_name = ['Fig ', datestr(now,'dd mmmm yyyy HH_MM_SS')];
savefig(Fig_name)


% Save Stage 2 results
% Re-run the DNA to FoM function to return the output state and heralding probability as well
[FigureOfMeritValue, output_density_operator, heralding_probability, FigureofMeritUnaugmentedValue, nbar_post_measurement] = DNAtoFigureOfMerit_s2_Anon(WinningDNA);
SaveResults_Paul_P2 (WinningDNA, FigureOfMeritValue, output_density_operator, heralding_probability, SelectedToolbox, num_modes, max_ops, Trunc_s2, output_loss, SelectedFigureOfMerit, use_heralding_probability, FigureofMeritUnaugmentedValue, nbar_post_measurement, exitflag, ga_output, final_population, final_scores, Cross_s2)




% Select the best DNA for stage 3
Corresponding_FOM_p2_results = final_scores;
Corresponding_FOM_ordered = sort(Corresponding_FOM_p2_results);    % Put the FOMs in order
Indices_to_keep = find(Corresponding_FOM_p2_results <= Corresponding_FOM_ordered(Population_size_s3),Population_size_s3);   % Get indices for the top DNA
DNA_for_stage3 = final_population(Indices_to_keep,:); % Select DNA to keep for next stage




%% Stage 3: Run the slower (but stricter? ga, or at least just smaller population) (maybe focussing on more local search)

% Load selected FOM for Bayes
load('FOM_Bayes_upto05.mat')

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

% Don't exceed time limit
ga_options.TimeLimit = Total_time_hrs*3600 - toc;

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
title('Stage 3')


% Re-run the DNA to FoM function to return the output state and heralding probability as well
[FigureOfMeritValue, output_density_operator, heralding_probability, FigureofMeritUnaugmentedValue, nbar_post_measurement] = DNAtoFigureOfMerit_s3_Anon(WinningDNA);
% Save results
SaveResults_Paul (WinningDNA, FigureOfMeritValue, output_density_operator, heralding_probability, SelectedToolbox, num_modes, max_ops, trunc, output_loss, SelectedFigureOfMerit, use_heralding_probability, FigureofMeritUnaugmentedValue, nbar_post_measurement, exitflag, ga_output, final_population, final_scores, Cross_s3)




%% Finalse

% Re-open the ga figure for the 1st ga, then delete the file
openfig(Fig_name);
Fig_full_name = [Fig_name,'.fig']; 
delete(Fig_full_name);

%
datetime('now')
toc


