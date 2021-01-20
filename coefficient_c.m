function c = coefficient_c(xs)
 ng=length(xs.D);
% DIFFUSION COEFFICIENT D    
    c=zeros(ng*6,1); % c (diffusion coefficient, D)  in 6N form : check Page2-112 of MATLAB PDETOOL for 3D Geometry) 
    % c = [D1 0 D1 0 0 D1 D2 0 D2 0 0 D2 D3 0 D3 0 0 D3... D_ng 0 D_ng 0 0 D_ng]
    for ig=1:1:ng  % For 2group diffusion, it wil stop at ng=2
        tmp = (ig-1)*6 +1;
        c(tmp)=xs.D(ig);  
        c(tmp+1)=0;
        c(tmp+2)=xs.D(ig); 
        c(tmp+3)=0;
        c(tmp+4)=0;
        c(tmp+5)=xs.D(ig);
    end
end