% clear the memory and the screen
clear; clc;

% =========================================================================
% Problem definition
% exact solution
exact = @(x) x.^3;

f = @(x) -6.0 * x;
g = 1;
h = 0;
% =========================================================================

% parameters of the FEM
n_el  = 13;       % number of elements
n_en  = 2;        % number of element nodes
n_int = 3;        % number of quadrature points
n_np  = n_el + 1; % number of points
n_eq  = n_np - 1; % number of equations

% =========================================================================
% Generate the mesh
% nodal coordinates
hh     = 1 / n_el;
x_coor = 0 : hh : 1;

% IEN
IEN = zeros(n_en, n_el);

for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(aa,ee) = ee + aa - 1;
  end
end
% =========================================================================

% ID
ID = 1 : n_np;
ID(end) = 0;

% LM
LM = ID(IEN);

% generate the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

K = spalloc(n_eq, n_eq, 3*n_eq);  % allocate the global stiffness matrix
F = zeros(n_eq, 1);     % allocate the global load vector

% Assembly of K and F
for ee = 1 : n_el
    k_e = zeros(n_en, n_en);
    f_e = zeros(n_en, 1);
    
    x_ele = x_coor(IEN(1:n_en,ee)); % A = IEN(a,e) and x_ele(a) = x_coor(A)

    k_e = zeros(n_en,n_en);
    f_e = zeros(n_en,1);

    x_ele = zeros(n_en,1);
    for a = 1 : n_en
        x_ele(a) = x_coor(IEN(a,e)); % A = IEN(a,e)
    end

    for l = 1 : n_int
        dx_dxi = 0.0;
        x_l = 0.0;
        for a = 1 : n_en
            dx_dxi = dx_dxi + x_ele(a) * PolyShape(a, xi(l), 1);
            x_l = x_l + x_ele(a) * PolyShape(a, xi(l), 0);
        end
        dxi_dx = 1.0 / dx_dxi;

        for a = 1 : n_en
            for b = 1 : n_en
                k_e(a,b) = k_e(a,b) + weight(l) * PolyShape(a, xi(l), 1) * PolyShape(b, xi(l), 1) * dxi_dx;
            end
        end

        for a = 1 : n_en
            f_e(a) = f_e(a) + weight(l) * PolyShape(a, xi(l), 0) * f(x_l) * dx_dxi;
        end

    end

    % Now we need to put element k and f into global K and F
    for a = 1 : n_en
        PP = LM(a,e);
        if PP > 0
            F(PP) = F(PP) + f_e(a);
            for b = 1 : 2
                QQ = LM(b,e);
                if QQ > 0
                    K(PP,QQ) = K(PP,QQ) + k_e(aa,bb);
                else
                    F(PP) = F(PP) - k_e(aa,bb) * g;
                end
            end
        end
    end

    % Neumann BC
    if e == 1
        F(ID(IEN(1,e))) = F(ID(IEN(1,e))) + h;
    end
end

% Now we have K and F assembled and we solve the linear system Kd = F
d_temp = K \ F;

disp = [uh; g];

% Check mid-point
e_mid = median(1:n_el);
Uh_x = disp(IEN(1, e_mid)) * PolyShape(1, 0, 1) * 2 / hh + disp(IEN(2, e_mid)) * PolyShape(2, 0.0, 1) * 2 / hh;
u_x = 3 * 0.5 * 0.5;

e_x = Uh_x - u_x
hh

% eof
