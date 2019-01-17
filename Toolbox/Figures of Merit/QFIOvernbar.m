function output = QFIOvernbar(density_operator)

    QFIVal = QFI(density_operator);
    nbar = Findnbar(density_operator,1);

    output = QFIVal/nbar;
    
    %% Avoid NaN by setting output=0 if nbar is zero
    if nbar < 1e-4
        output=0;
    end
    
end