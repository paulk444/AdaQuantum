function [state] = SingleModeSqueezedState(z, trunc)
% Creates a squeezed state. Takes [|z|, arg(z)] as argument

r = z(1);
phi = z(2);

state = sparse(trunc,1);

for n = 0:floor((trunc-1)/2)
    state = state + sqrt(factorial(2*n))/factorial(n) * (- 1/2)^n * exp(1i * n * phi) * tanh(r)^n * FockState(2*n, trunc);
end


state = sqrt(sech(r)) * state;

end