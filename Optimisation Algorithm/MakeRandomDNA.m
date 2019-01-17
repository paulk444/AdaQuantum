function DNA = MakeRandomDNA(LowerBounds, UpperBounds, Integers)

% Makes a random DNA strand using the LowerBounds, UpperBounds, and whether or not each DNA element is an integer or not

rng('shuffle')  % Make the random numbers different each time

    DNA_length = length(LowerBounds);
    
    DNA = rand(DNA_length, 1);
    
    for i = 1:DNA_length
       if ismember(i, Integers)
           DNA(i) = IntegerPicker(DNA(i), LowerBounds(i), UpperBounds(i));
       else
           DNA(i) = (UpperBounds(i) - LowerBounds(i)) * DNA(i) + LowerBounds(i);
       end
    end

end