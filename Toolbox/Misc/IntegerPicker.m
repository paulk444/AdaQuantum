function SelectedInt = IntegerPicker (selector, min_int, max_int)
% Takes a selector (value between 0 and 1) and uses it to pick an integer
% between min_int and max_int
% Note: need 0<selector<1

    SelectedInt = min_int + round((max_int - min_int + 1) * selector - 0.5);


end