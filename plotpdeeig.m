function plotpdeeig(model,result)
    
    u = result.Eigenvectors;

    if sum(sum(u(:,1,1))) <0
        u=-u; 
    end

    figure
    pdeplot3D(model,'ColorMap',u(:,1,1))
    daspect([max(daspect)*[1 1] 2]) 
    grid on
    title(['k =' num2str(1/result.Eigenvalues(1),6) , '. Group 1 flux.'])
    
    figure
    pdeplot3D(model,'ColorMap',u(:,2,1))
    daspect([max(daspect)*[1 1] 2]) 
    grid on
    title(['k =' num2str(1/result.Eigenvalues(1),6) , '. Group 2 flux.'])
    disp(['k = ' num2str(1/result.Eigenvalues(1),8) '    Eigenvalue = ' num2str(result.Eigenvalues(1))]);
end