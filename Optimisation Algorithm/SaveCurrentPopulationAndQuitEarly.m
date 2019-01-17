function [state,options,optchanged] = SaveCurrentPopulationAndQuitEarly(options,state,~)
% ONLY VALID FOR T2 SETTINGS    


%every 10 generations, save the population to a variable in the base workspace
    if rem(state.Generation,10) == 0
        
        current_population = state.Population;
        assignin('base','last_saved_population', current_population);
        assignin('base','last_saved_generation', state.Generation);
        
        % If both of the operation selctors have below 10% displacement operators, quit
        if length(find(current_population(:,12)==2))/length((current_population(:,12))) < 0.13 && length(find(current_population(:,7)==2))/length((current_population(:,7))) < 0.13
            state.StopFlag = 'y';
        end
        
    end
    
    optchanged = false;    
    
end