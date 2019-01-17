function [Nmat] = NumberOperator (trunc, acting_on_mode, number_of_modes)

    Nmat = CreationOperatorGeneral(trunc, acting_on_mode, number_of_modes) * AnnihilationOperatorGeneral(trunc, acting_on_mode, number_of_modes);

end