function val = PolyShape(a, xi, der)

if a == 1
    if der == 0
        val = 0.5 * (1-xi);
    elseif der == 1
        val = -0.5;
    end
elseif a == 2
    if der == 0
        val = 0.5 * (1+xi);
    elseif der == 1
        val = 0.5;
    else
    end
end