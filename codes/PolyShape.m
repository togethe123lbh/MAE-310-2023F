function val = PolyShape(a, xi, der)

%linear element
% if a == 1
%     if der == 0
%         val = 0.5 * (1-xi);
%     elseif der == 1
%         val = -0.5;
%     end
% elseif a == 2
%     if der == 0
%         val = 0.5 * (1+xi);
%     elseif der == 1
%         val = 0.5;
%     else
%     end
% end

% quadratic element
if a==1
    if der == 0
        val=0.5*xi*(xi-1);
    elseif der == 1
        val=-0.5+xi;
    end
elseif a == 2
    if der==0
        val=1-xi^2;
    elseif der==1
        val=-2*xi;
    end
elseif a == 3
    if der==0
        val=0.5*xi*(1+xi);
    elseif der==1
        val=0.5+xi;
    end
end

%cubic element
% if a == 1
%     if der == 0
%         val = 0.5 * (1-xi);
%     elseif der == 1
%         val = -0.5;
%     end
% elseif a == 2
%     if der == 0
%         val = 0.5 * (1+xi);
%     elseif der == 1
%         val = 0.5;
%     else
%     end
% elseif a == 3
%     if der == 0
%         val = 0.5 * (1-xi);
%     elseif der == 1
%         val = -0.5;
%     end
% elseif a == 4
%     if der == 0
%         val = 0.5 * (1+xi);
%     elseif der == 1
%         val = 0.5;
%     else
%     end
% end


