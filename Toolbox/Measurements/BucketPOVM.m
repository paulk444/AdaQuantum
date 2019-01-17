function POVM = BucketPOVM(detect_photons, detector_loss, trunc)

% Returns a POVM for a Bucket detector with loss 'detector_loss'
% Options for detect_photons:   false/0 measure |0><0| (i.e. no photons)
%                               true/1 measure I - |0><0| (i.e. at least one photon)

%===== Lana Mineh, Rosanna Nichols, and  Paul Knott 2018 =========%
%============= https://arxiv.org/abs/1812.01032 ==================%

    % Construct the POVM E_0
    E_0 = 0;
    for n = 0:trunc-1
        Fock = FockState(n,trunc);
        E_0 = E_0 + detector_loss^n*(Fock*Fock');
    end
    
    if detect_photons
        POVM = speye(trunc) - E_0;
    else
        POVM = E_0;
    end
        
end