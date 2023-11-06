% clear the memory and the screen
clear all; clc;

% exact solution
exact = @(x) x.^3;

% problem definition
f = @(x) -6.0 * x;
g = 1;
h = 0;

% generate my mesh
n_el = 9;
hh = 1 / n_el;
x_coor = 0 : hh : 1;

% number of element nodes
n_en = 2;

% IEN
IEN = zeros(n_en, n_el);

for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(aa,ee) = ee + aa - 1;
  end
end

% ID
n_pt = n_el + 1; % number of points
ID = 1 : n_pt;
ID(end) = 0;

% LM
LM = ID(IEN);

n_eq = n_pt - 1; % number of equations

% generate the quadrature rule
n_int = 30;
[xi, weight] = Gauss(n_int, -1, 1);

K = spalloc(n_eq, n_eq, 3*n_eq); % allocate the global stiffness matrix
F = zeros(n_eq, 1);    % allocate the global load vector

% Assembly of K and F
for ee = 1 : n_el

    k_e = zeros(n_en, n_en);
    f_e = zeros(n_en, 1);

    x_ele = zeros(n_en,1);
    for aa = 1 : n_en
        x_ele(aa) = x_coor(IEN(aa,ee)); % A = IEN(a,e)
    end

    for l = 1 : n_int
        dx_dxi = 0.0;
        x_l = 0.0;
        for aa = 1 : n_en
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(aa, xi(l), 1);
            x_l = x_l + x_ele(aa) * PolyShape(aa, xi(l), 0);
        end
        dxi_dx = 1.0 / dx_dxi;

        for aa = 1 : n_en
            for bb = 1 : n_en
                k_e(aa,bb) = k_e(aa,bb) + weight(l) * PolyShape(aa, xi(l), 1) * PolyShape(bb, xi(l), 1) * dxi_dx;
            end
        end

        for aa = 1 : n_en
            f_e(aa) = f_e(aa) + weight(l) * PolyShape(aa, xi(l), 0) * f(x_l) * dx_dxi;
        end

    end

    % Now we need to put element k and f into global K and F
    for aa = 1 : n_en
        PP = LM(aa,ee);
        if PP > 0
            F(PP) = F(PP) + f_e(aa);
            for bb = 1 : n_en
                QQ = LM(bb,ee);
                if QQ > 0
                    K(PP,QQ) = K(PP,QQ) + k_e(aa,bb);
                else
                    F(PP) = F(PP) - k_e(aa,bb) * g;
                end
            end
        end
    end

    if ee == 1
        F(ID(IEN(1,ee))) = F(ID(IEN(1,ee))) + h;
    end
end

% Now we have K and F
% Solve Kd = F
uh = K \ F;

disp = [uh; g];

% Post processing-1 error analysis of the derivative error at x = 0.5
exact_x = @(x) 3.0 * x.^2;
e_m = median(1:n_el);
A_m = IEN(1, e_m);
uh_x = disp(A_m) * PolyShape(1, 0.0, 1) * (2 / hh) ...
  + disp(A_m+1) * PolyShape(2, 0.0, 1) * (2 / hh);

e_x = uh_x - exact_x(0.5);




% eof