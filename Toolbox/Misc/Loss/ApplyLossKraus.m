function output_state = ApplyLossKraus(input_state, loss_rate, cut_off)
% Function to apply loss using Kraus operators
% Takes density matrix as input
% Cut_off (cut_off \leq trunc) defines how many Kraus operators to apply - how many photons likely to be lost

%=========== Lana Mineh and Paul Knott 2018 ================%
%========== https://arxiv.org/abs/1812.01032 ===============%

trunc = length(input_state);
output_state = sparse(trunc,trunc);

for k = 0:cut_off-1
    E_k = KrausOperator(loss_rate, k, trunc);
    output_state = output_state + E_k*input_state*E_k'; 
end

% If trunc is too large then the factorial gives Inf
    if trunc > 171
        disp('WARNING: Truncation MIGHT BE too large for ApplyLossKraus, as factorial(172) gives Inf')
    end
    
end

%% Function to calculate the Kraus Operators

function E_k = KrausOperator(loss_rate, k, trunc)

E_k = sparse(trunc, trunc);

for n = k:trunc-1
    E_k(n-k+1,n+1) = sqrt(factorial(n))*1/sqrt(factorial(k))*1/sqrt(factorial(n-k))*sqrt((1-loss_rate)^(n-k)*loss_rate^k);
end
    
end
