function d = coefficient_d(xs)
ng=length(xs.vf);
% miu-fission and Chi
    d=zeros(ng,ng);% For 2group diffusion, 'd' will be 2 x 2 zero matrix
    for ig=1:1:ng
        d(ig,:) = xs.X(ig)*xs.vf(:);  % For 2group role 1; ??1f ??2f. Note: xs.X1=1 is fission fraction & X2=0
        continue
    end
    d=d(:);
end
