function SaveResults (WinningDNA, FigureOfMeritValue, output_density_operator, heralding_probability, SelectedToolbox, num_modes, max_ops, trunc, output_loss, SelectedFigureOfMerit, use_heralding_probability, FigureofMeritUnaugmentedValue, nbar_post_measurement, exitflag, solver_output, final_population, final_scores, crossover_fraction, solver_options,use_algorithm)

% make new folder based on current date and time
if FigureofMeritUnaugmentedValue > 30
    FolderName = ['good one ', use_algorithm, 'results ', datestr(now,'dd mmmm yyyy HH_MM')];
else
    FolderName = [use_algorithm, ' results ', datestr(now,'dd mmmm yyyy HH_MM')];
end
disp(FolderName)
mkdir(FolderName)
% save MATLAB variables
save([FolderName, '/DNA_save'], 'WinningDNA', 'FigureOfMeritValue', 'output_density_operator', 'heralding_probability', 'SelectedToolbox', 'num_modes', 'max_ops', 'trunc', 'output_loss', 'SelectedFigureOfMerit', 'use_heralding_probability', 'FigureofMeritUnaugmentedValue', 'nbar_post_measurement', 'exitflag', 'solver_output', 'final_population', 'final_scores', 'crossover_fraction','solver_options')
% save genetic algorithm figure
savefig([FolderName, '/solver plots.fig'])
% print DNA to command window whilst recording the console in a text file
diary([FolderName, '/DNA print.txt'])
DNAtoDisplay(WinningDNA, SelectedToolbox, SelectedFigureOfMerit, num_modes, max_ops);
disp(FigureOfMeritValue)
diary off

end