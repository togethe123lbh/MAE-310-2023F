% clear the memory and the screen
clear; clc;

% =========================================================================
% Problem definition exact solution

%linear and quadratic
% g = 1;
% h = 0;
% exact = @(x) x.^3;
% exact_dx=@(x) 3*(x.^2);
% f = @(x) -6.0 * x;
% exact_2=@(x) x.^6;
% el2down=integral(exact_2,0,1);
% exact_dx2=@(x) 9*(x.^4);
% eH2down=integral(exact_dx2,0,1);

%cubic
g = sin(1);
h = -1;
exact = @(x) sin(x);
exact_dx=@(x) cos(x);
handle = @(x) sin(x);
exact_2=@(x) sin(x).^2;
el2down=integral(exact_2,0,1);
exact_dx2=@(x) cos(x).^2;
eH2down=integral(exact_dx2,0,1);
% =========================================================================

% parameters of the FEM
n_el  = 2;       % number of elements
%linear
%n_en  = 2;        % number of element nodes

%quadratic
% n_en=3;

%cubic
n_en=4;
      % number of quadrature points

%for el=2:2:16
 el  = n_el;
n_np  = 3*el + 1; % number of points 
n_eq  = n_np - 1; % number of equations

% =========================================================================
% Generate the mesh nodal coordinates
hh     = 1 / (3*el);
x_coor = 0 : hh : 1;

% IEN
IEN = zeros(n_en, n_el);

for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(aa,ee) = 3*ee + aa - 3;
  end
end
% =========================================================================

% ID and LM arrays are generated based on the BC info
ID = 1 : n_np;
ID(end) = 0; % Modify ID according to the Dirichlet BC info

LM = ID(IEN);
el2=zeros(6,1);
eH2=zeros(6,1);
% generate the quadrature rule
% for n_int=1:1:6

% n_int = 1;
Final=zeros(6,1);
for n_int=1:6
[xi, weight] = Gauss(n_int, -1, 1);

K = spalloc(n_eq, n_eq, 3*n_eq); % allocate the global stiffness matrix
K_exact=spalloc(n_eq, n_eq, 3*n_eq);
F = zeros(n_eq, 1);    % allocate the global load vector
F_exact= zeros(n_eq, 1); 
% Assembly of K and F
for ee = 1 : n_el
    k_e = zeros(n_en, n_en);
    k_exact=zeros(n_en,n_en);
    f_e = zeros(n_en, 1); 
    f_exact = zeros(n_en, 1);
    x_ele = x_coor(IEN(1:n_en,ee)); % A = IEN(a,e) and x_ele(a) = x_coor(A)
dx_dxiexact=@(x) 0;
    for ll = 1 : n_int
        dx_dxi = 0.0;
        x_l = 0.0;

        for aa = 1 : n_en
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(aa, xi(ll), 1);
            
            dx_dxiexact=@(x) dx_dxiexact(x)+ x_ele(aa) *PolyShape(aa, x, 1);
            x_l = x_l + x_ele(aa) * PolyShape(aa, xi(ll), 0);
        end
        dxi_dx = 1.0 /dx_dxi;
      
       

        for aa = 1 : n_en
            f_e(aa) = f_e(aa) + weight(ll) * PolyShape(aa, xi(ll), 0) * handle(x_l) * dx_dxi;
            f_exact1=@(x) PolyShape(aa,x,0).*handle(x_l).* dx_dxi;
            f_exact(aa)=integral(f_exact1,-1,1);
            for bb = 1 : n_en
                k_e(aa,bb) = k_e(aa,bb) + weight(ll) * PolyShape(aa, xi(ll), 1) * PolyShape(bb, xi(ll), 1) * dxi_dx;
              %  k_exact(aa,bb)=integral(PolyShape(aa, x, 1) * PolyShape(bb, x, 1)*dxi_dx,-1,1);
              handle=@(x) PolyShape(aa, x, 1) .* PolyShape(bb, x, 1).*dxi_dx;
              k_exact(aa,bb)=integral(handle,-1,1);
              
            end
        end

%         for bb = 1 : n_en
%             for aa = 1 : n_en
%                 for cc = 1 : n_en
%             dxi_dx=dxi_dx+PolyShape(cc, x, 1);
%                 end
%                 f=@(x) PolyShape(aa, x, 1) .* PolyShape(bb, x, 1);
% k_exact(aa,bb)=k_exact(aa,bb)+integral(f,-1,1);
%             end
%         end
    end

    % Now we need to put element k and f into global K and F
    for aa = 1 : n_en
        PP = LM(aa,ee);
        if PP > 0
            F(PP) = F(PP) + f_e(aa);
            F_exact(PP)=F_exact(PP) +f_exact(aa);
            for bb = 1 : n_en
                QQ = LM(bb,ee);
                if QQ > 0
                    K(PP,QQ) = K(PP,QQ) + k_e(aa,bb);
                    K_exact(PP,QQ)=K_exact(PP,QQ)+k_exact(aa,bb);
                else
                    F(PP) = F(PP) - k_e(aa,bb) * g;
                     F_exact(PP) = F_exact(PP) - k_exact(aa,bb) * g;
                end
            end
        end
    end

    if ee == 1
        F(ID(IEN(1,ee))) = F(ID(IEN(1,ee))) + h;
        F_exact(ID(IEN(1,ee))) = F_exact(ID(IEN(1,ee))) + h;
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
exact_M=0.0;
exact_dxM=0.0;
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
%         el2up1=el2up1+weight(ll)*(uh-exact(xl)).^2*dx_dxi;
%        % el2down1=el2down1+weight(ll)*(exact(xl)).^2*dx_dxi;
%         eH2up1=eH2up1+weight(ll)*(uh_dx*dxi_dx-exact_dx(xl)).^2*dx_dxi;
%        % eH2down1=eH2down1+(exact_dx(xl)).^2*dx_dxi*weight(ll);
% exact_M=exact_M+weight(ll)*(exact(xl)).^2*dx_dxi;
% exact_dxM=exact_dxM+(exact_dx(xl)).^2*dx_dxi*weight(ll);

    end
end
fi=0;
for aa=1:n_eq
fi=fi+(F(aa)-F_exact(aa))*(F(aa)-F_exact(aa));
end

Final(n_int)=sqrt(fi);
end

% for aa=3:20
% plot(aa,Final(aa),'o');
% hold on
% end
% eH2(n_int,1)=sqrt(exact_dxM-eH2down);
% % end
%  
% x=1:1:6;
% plot(x,el2)
% hold on
% plot(x,eH2)
% ele_hh=ones(8,6)
% 1./(2:2:16);
% for col=1:1:6
% for row=1:1:8
%  plot(log(ele_hh),log(el2(row,col)))
%  hold on
%  plot(log(ele_hh),log(eH2(row,col)))
%  slope_el2=(log(el2(n_el*0.5,end))-log(el2(n_el*0.5,1)))/(log(ele_hh(end))-log(ele_hh(1)));
%  slope_eH2=(log(eH2(n_el*0.5,end))-log(eH2(n_el*0.5,1)))/(log(ele_hh(end))-log(ele_hh(1)));
% end   
%  end