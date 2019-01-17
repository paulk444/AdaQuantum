function output = TensorProduct(things)

output = 1;

for i = 1:length(things)
    
    output = kron(output, things{i});
    
end

end