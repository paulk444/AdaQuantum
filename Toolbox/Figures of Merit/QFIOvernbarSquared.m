function output = QFIOvernbarSquared(density_operator)

    QFIVal = QFI(density_operator);
    nbar = Findnbar(density_operator,1);

    output = QFIVal/(nbar^2);
    
    %% Avoid NaN by setting output=0 if nbar is zero
    if nbar < 1e-4
        output=0;
    end
    
end