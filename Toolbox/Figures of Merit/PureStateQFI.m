function fisherinfo = PureStateQFI(density_operator)

trunc = length(density_operator);
N = NumberOperator(trunc, 1, 1);
NSquared = N^2;

fisherinfo = 4 * (trace(NSquared * density_operator) - abs(trace(N * density_operator))^2);


end