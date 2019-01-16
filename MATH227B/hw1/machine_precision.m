function epsilon = machine_precision()

    epsilon = 1.0;
    while 1.0 + (epsilon/2) ~= 1.0
        epsilon = (epsilon / 2);
        
end