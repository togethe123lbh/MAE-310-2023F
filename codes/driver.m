% clear the memory and the screen
clear; clc;

% =========================================================================
% Problem definition exact solution


g = 1;
h = 0;
exact = @(x) x.^3;
exact_dx=@(x) 3*(x.^2);
f = @(x) -6.0 * x;
exact_2=@(x) x.^6;
el2down=integral(exact_2,0,1);
exact_dx2=@(x) 9*(x.^4);
eH2down=integral(exact_dx2,0,1);
% =========================================================================

% parameters of the FEM
n_el  = 2;       % number of elements
n_en  = 2;        % number of element nodes
n_int = 3;        % number of quadrature points
el2=zeros(8,1);
eH2=zeros(8,1);
for el=2:2:16
    n_el = el;
n_np  = n_el + 1; % number of points 
n_eq  = n_np - 1; % number of equations

% =========================================================================
% Generate the mesh nodal coordinates
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

% ID and LM arrays are generated based on the BC info
ID = 1 : n_np;
ID(end) = 0; % Modify ID according to the Dirichlet BC info

LM = ID(IEN);

% generate the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

K = spalloc(n_eq, n_eq, 3*n_eq); % allocate the global stiffness matrix
F = zeros(n_eq, 1);    % allocate the global load vector

% Assembly of K and F
for ee = 1 : n_el
    k_e = zeros(n_en, n_en);
    f_e = zeros(n_en, 1); 
    
    x_ele = x_coor(IEN(1:n_en,ee)); % A = IEN(a,e) and x_ele(a) = x_coor(A)

    for ll = 1 : n_int
        dx_dxi = 0.0;
        x_l = 0.0;
        for aa = 1 : n_en
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(aa, xi(ll), 1);
            x_l = x_l + x_ele(aa) * PolyShape(aa, xi(ll), 0);
        end
        dxi_dx = 1.0 / dx_dxi;

        for aa = 1 : n_en
            f_e(aa) = f_e(aa) + weight(ll) * PolyShape(aa, xi(ll), 0) * f(x_l) * dx_dxi;
            for bb = 1 : n_en
                k_e(aa,bb) = k_e(aa,bb) + weight(ll) * PolyShape(aa, xi(ll), 1) * PolyShape(bb, xi(ll), 1) * dxi_dx;
            end
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

% Now we have K and F assembled and we solve the linear system Kd = F
d_temp = K \ F;

% Generate the full solution vector by inserting back the Dirichlet value
disp = [d_temp; g];

% eof
%2b
el2up1=0.0;
el2down1=0.0;
eH2up1=0.0;
eH2down1=0.0;
for ee=1:n_el
    for ll=1:n_int
        uh=0.0;
        uh_dx=0.0;
        xl=0.0;
        dx_dxi=0.0;
        for aa=1:n_en
            uh=uh+disp(IEN(aa,ee))*PolyShape(aa,xi(ll),0);
            uh_dx=uh_dx+disp(IEN(aa,ee))*PolyShape(aa,xi(ll),1);
            xl=xl+x_coor(IEN(aa,ee))*PolyShape(aa,xi(ll),0);
            dx_dxi=dx_dxi+x_coor(IEN(aa,ee))*PolyShape(aa,xi(ll),1);
        end
        dxi_dx= 1.0 / dx_dxi;
        el2up1=el2up1+weight(ll)*(uh-exact(xl)).^2*dx_dxi;
       % el2down1=el2down1+weight(ll)*(exact(xl)).^2*dx_dxi;
        eH2up1=eH2up1+weight(ll)*(uh_dx*dxi_dx-exact_dx(xl)).^2*dx_dxi;
       % eH2down1=eH2down1+(exact_dx(xl)).^2*dx_dxi*weight(ll);
    end
end
el2(n_el*0.5)=sqrt(el2up1/el2down);
eH2(n_el*0.5)=sqrt(eH2up1/eH2down);
end

    ele_hh=1./(2:2:16);
 plot(log(ele_hh),log(el2))
 hold on
 plot(log(ele_hh),log(eH2))
 slope_el2=(log(el2(end))-log(el2(1)))/(log(ele_hh(end))-log(ele_hh(1)));
 slope_eH2=(log(eH2(end))-log(eH2(1)))/(log(ele_hh(end))-log(ele_hh(1)));
