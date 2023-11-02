% clear the memory and the screen
clear all; clc;

% exact solution
exact = @(x) sin(x);

% problem definition
f = @(x) sin(x);
g = sin(1);
h = -cos(0);

% generate my mesh
n_el = 5;
hh = 1 / n_el;
x_coor = 0 : hh : 1;

% IEN
IEN = zeros(2, n_el);

for e = 1 : n_el
    IEN(1,e) = e;
    IEN(2,e) = e+1;
end

% ID
n_pt = n_el + 1; % number of points
ID = 1 : n_pt;
ID(end) = 0;

n_eq = n_pt - 1; % number of equations

% generate the quadrature rule
n_int = 30;
[xi, weight] = Gauss(n_int, -1, 1);

K = zeros(n_eq, n_eq); % allocate the global stiffness matrix
F = zeros(n_eq, 1);    % allocate the global load vector

% Assembly of K and F
for e = 1 : n_el

    k_e = zeros(2,2);
    f_e = zeros(2,1);

    x_ele = zeros(2,1);
    for a = 1 : 2
        x_ele(a) = x_coor(IEN(a,e)); % A = IEN(a,e)
    end

    for l = 1 : n_int
        dx_dxi = 0.0;
        x_l = 0.0;
        for a = 1 : 2
            dx_dxi = dx_dxi + x_ele(a) * PolyShape(a, xi(l), 1);
            x_l = x_l + x_ele(a) * PolyShape(a, xi(l), 0);
        end
        dxi_dx = 1.0 / dx_dxi;

        for a = 1 : 2
            for b = 1 : 2
                k_e(a,b) = k_e(a,b) + weight(l) * PolyShape(a, xi(l), 1) * PolyShape(b, xi(l), 1) * dxi_dx;
            end
        end

        for a = 1 : 2
            f_e(a) = f_e(a) + weight(l) * PolyShape(a, xi(l), 0) * f(x_l) * dx_dxi;
        end

    end

    % Now we need to put element k and f into global K and F
    for a = 1 : 2
        AA = IEN(a,e);
        PP = ID(AA);
        if PP > 0
            F(PP) = F(PP) + f_e(a);
            for b = 1 : 2
                BB = IEN(b,e);
                QQ = ID(BB);
                if QQ > 0
                    K(PP,QQ) = K(PP,QQ) + k_e(a,b);
                else
                    F(PP) = F(PP) - k_e(a,b) * g;
                end
            end
        end
    end

    if e == 1
        F(ID(IEN(1,e))) = F(ID(IEN(1,e))) + h;
    end
end

% Now we have K and F
% Solve Kd = F
uh = K \ F;

d = [uh; g];

% eof