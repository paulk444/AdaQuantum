function POVM = MultiplexDetectorPOVM(no_to_detect, detector_settings, trunc)
% Returns a POVM for a number measurement using a multiplex detector, detecting 'no_to_detect' photons
% detector settings = [detectors, detector_loss]
    % detectors = number of detectors 
    % no_to_detect = the POVM element E_(no_to_detect) you want to measure
    
%=========== Lana Mineh and Paul Knott 2018 ================%
%========== https://arxiv.org/abs/1812.01032 ===============%


    detectors = detector_settings(1);
    detector_loss = detector_settings(2);
    
    POVM = 0;
    
    % Construct POVM
    for c = 0:trunc-1
        if c >= no_to_detect
            w = factorial(detectors)*Stirling2(c,no_to_detect)/(factorial(detectors-no_to_detect)*detectors^c);
            E_c = NumberPOVM(c, detector_loss, trunc);
            POVM = POVM + w*E_c;
        end
    end

end

%% Calculate Stirling number of the second kind

function S = Stirling2(n, k)

S = 0;
for j = 0:k
    S = S + (-1)^(k-j)*1/factorial(j)*1/factorial(k-j)*j^n;
end

end

