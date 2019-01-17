function fisherinfo = QFI(density_operator)

isPure = PurityTest(density_operator, 0.001);

if isPure
    fisherinfo = PureStateQFI(density_operator);
else
    fisherinfo = MixedStateQFI(density_operator);
end
% Note: both pure and mixed QFI take 1-mode state, apply exp(i nhat phi),
% then calculate QFI. i.e. they don't consider phase difference in an
% interferometer.


end