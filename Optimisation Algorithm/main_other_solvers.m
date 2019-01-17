%% Load toolbox + UI setting selection

% shuffle random numbers, to be safe
rng('shuffle')

Toolbox = 'Regular'; % options: 'OAM' or 'Regular' (our full quantum optics one)
use_UI = true;
use_integers = false;
use_random_algorithm = false;
use_algorithm_options = {'ga', 'swarm', 'pattern_search', 'simulated_annealing'};
use_algorithm = use_algorithm_options{4}; % 1, 2, 3, 4 give 'ga', 'swarm', 'pattern_search', 'simulated_annealing', respectively

heralding_accepted = 1e-3;

if strcmp( Toolbox, 'Regular')
    %load default toolbox
    load('Full Toolbox.mat')
elseif strcmp( Toolbox, 'OAM')
    %load OAM toolbox
    load('OAM Toolbox.mat')
else
    disp('Error: no toolbox loaded')
end

if use_UI
    
    % make UI for toolbox selection, return user selections in selected toolbox
    [SelectedToolbox, num_modes, max_ops, trunc, output_loss]  = MakeUI (FullToolbox, ToolboxTypes, default_num_modes, default_max_ops, default_trunc, default_output_loss);
    % display the settings selected
    MakeSelectedUI (SelectedToolbox, ToolboxTypes, num_modes, max_ops, trunc, output_loss);
    % make UI for figure of merit selection, return user selections
    [SelectedFigureOfMerit, use_heralding_probability] = MakeFigureOfMeritUI (FiguresOfMerit);
    % display the FoM settings selected
    MakeSelectedFigureOfMeritUI (SelectedFigureOfMerit, use_heralding_probability)
    
else
    if strcmp( Toolbox, 'Regular')
        % load T2 settings
        load('T2_settings.mat')
        load('5594settings.mat')
    else
        SelectedToolbox = FullToolbox;
        num_modes = default_num_modes;
        max_ops = default_max_ops;
        trunc = 3;
        output_loss = default_output_loss;
        SelectedFigureOfMerit = FiguresOfMerit;
        use_heralding_probability = false;
        SelectedFigureOfMerit{2}{5} = feval(SelectedFigureOfMerit{2}{5}, SelectedFigureOfMerit{2}{4});
    end
    
end

%% automated setting finding + adjusting

% Find the settings for the genetic algorithm given the toolbox selected from the UI.
[DNA_length, LowerBounds, UpperBounds, Integers] = FindDNASettings(SelectedToolbox, SelectedFigureOfMerit, num_modes, max_ops, use_integers);

% make and anonymous function to pass in the parameters that are not the DNA (ga allows only the DNA as an argument)
DNAtoFigureOfMeritAnon = @(DNA) DNAtoFigureOfMerit_newSafety(DNA, SelectedToolbox, SelectedFigureOfMerit, use_heralding_probability, num_modes, max_ops, trunc, output_loss, heralding_accepted);        
    
    if use_random_algorithm
        int_picker = rand;
        algorithm_picker = IntegerPicker(int_picker, 1, 4);
        use_algorithm = use_algorithm_options{algorithm_picker};
    end
    
    if strcmp( use_algorithm, 'ga' )
        
        solver_options = ga_options;
        
        % fiddle with options
        solver_options.PopulationSize = 100;
        solver_options.EliteCount = 10;
        solver_options.StallGenLimit = 50;
        %     ga_options.CrossoverFcn = @crossovertwopoint;
        solver_options.HybridFcn = {@patternsearch, patternsearch_options};
        
        if strcmp( Toolbox, 'Regular')
            
           solver_options.OutputFcns = {@SaveCurrentPopulationAndQuitEarly}; 
            
        end
        
        
    elseif strcmp( use_algorithm, 'swarm' )
        
        solver_options = swarm_options;
        
        solver_options.HybridFcn = {@patternsearch, patternsearch_options};
        solver_options.SwarmSize = 100;
        
    elseif strcmp( use_algorithm, 'pattern_search' )
        
        solver_options = patternsearch_options;
        
    elseif strcmp( use_algorithm, 'simulated_annealing' )
        
        solver_options = anneal_options;
        solver_options.TemperatureFcn = @temperatureboltz;
%         solver_options.HybridFcn = {@patternsearch, patternsearch_options};
        solver_options.PlotFcns = { @saplotbestf, @saplottemperature, @saplotf, @saplotstopping };
        solver_options.StallIterLimit = Inf;
        solver_options.MaxFunEvals = Inf;
        solver_options.MaxIter = Inf;
        solver_options.ReannealInterval = Inf;
        solver_options.InitialTemperature = 100;
        
    end
    
 while true
    
    %% run algorithm
    
    if strcmp(use_algorithm, 'ga')
        
        % set up for parallel processing (gcp is "get current pool" - if one is not already running it creates one)
        p = gcp;
        
        crossover_fraction = rand;
        ga_options.CrossoverFraction = crossover_fraction;
        
        % run genetic algorithm
        [WinningDNA, ~, exitflag, solver_output, final_population, final_scores] = ga(DNAtoFigureOfMeritAnon,DNA_length,[],[],[],[],LowerBounds,UpperBounds,[],Integers,solver_options);
        
    elseif strcmp( use_algorithm, 'swarm' )
        
        % set up for parallel processing (gcp is "get current pool" - if one is not already running it creates one)
        p = gcp;
        
        [WinningDNA, ~, exitflag, solver_output] = particleswarm(DNAtoFigureOfMeritAnon, DNA_length, LowerBounds, UpperBounds, solver_options);
        
    elseif strcmp( use_algorithm, 'pattern_search' )
        
        % set up for parallel processing (gcp is "get current pool" - if one is not already running it creates one)
        p = gcp;
        
        initial_point = MakeRandomDNA(LowerBounds, UpperBounds, []);
        [WinningDNA, ~, exitflag, solver_output] = patternsearch(DNAtoFigureOfMeritAnon, initial_point, [],[],[],[],LowerBounds,UpperBounds,[], solver_options);
        
    elseif strcmp( use_algorithm, 'simulated_annealing' )
        
        initial_point = MakeRandomDNA(LowerBounds, UpperBounds, []);
        [WinningDNA, ~, exitflag, solver_output] = simulannealbnd(DNAtoFigureOfMeritAnon, initial_point, LowerBounds, UpperBounds, solver_options);
        
    end
    
    if strcmp(use_algorithm, 'ga')
        % adjust figure axes to something sensible
        fig = gcf;
        fig.Children (5).XLim = [1 inf];
        fig.Children (7).XLim = [1 inf];
    end
    %
    
FolderName = ['Results, ', datestr(now,'dd mmmm yyyy HH_MM')];
mkdir(FolderName)
% Save ENTIRE workspace!
save([FolderName, '/Workspace'])
    
    
end