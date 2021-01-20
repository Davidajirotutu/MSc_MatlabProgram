function a = coefficient_a(xs)
ng=length(xs.a);
% Absorption (LHS of the equation) and scattering source (RHS of the equation)
    a=zeros(ng,ng); % For 2group diffusion, 'a' will be 2 x 2 zero matrix
    for ig=1:1:ng
        a(ig,ig) = xs.a(ig); 
        for igg=1:1:ng
           if igg~=ig  
               a(ig,igg) = -xs.s(igg,ig); % For 2group diffusion, a(2,1) will be -(?1->2) 
               a(ig,ig) = a(ig,ig) + xs.s(ig,igg);  % For 2group diffusion, 'a(1,1) will be ?a1 +(?1->2)
           end
        end
    end
    a=a(:); % change it to a column vector
end