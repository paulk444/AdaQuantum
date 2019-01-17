function isequal = equaltol (thing1, thing2, tol)

isequal = norm(thing1 - thing2) < tol;

if isequal
    disp('Correct')
else
    disp('ERROR')
end

end