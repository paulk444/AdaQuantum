function [state,options,optchanged] = SaveCurrentPopulation(options,state,~)
    
%every 10 generations, save the population to a variable in the base workspace
    if rem(state.Generation,10) == 0
        
        current_population = state.Population;
        assignin('base','last_saved_population', current_population);
        assignin('base','last_saved_generation', state.Generation);
        
    end
    
    optchanged = false;    
    
end