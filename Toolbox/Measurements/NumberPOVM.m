function POVM = NumberPOVM(number, detector_loss, trunc)

% Returns a POVM for a number measurement of 'number' photons, for a detector with loss 'detector_loss'
% number = number of photons heralded/post-selected

%=========== Lana Mineh and Paul Knott 2018 ================%
%========== https://arxiv.org/abs/1812.01032 ===============%

    % Construct POVM
    POVM = 0;
    for k = number:trunc-1
        Fock = FockState(k,trunc);
        POVM = POVM + factorial(k)/(factorial(number)*factorial(k-number))*detector_loss^(k-number)*(Fock*Fock');
    end
    POVM = (1-detector_loss)^number*POVM;
    
    % If trunc is too large then the factorial gives Inf
    if trunc > 171
        disp('WARNING: Truncation too large for NumberPOVM, as factorial(172) gives Inf')
    end
    
end