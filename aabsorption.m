function a = aabsorption(xs,xcoordinates,ycoordinates,rectarray,location)
np = numel(location.x);
n1 = 4;
zcoordinates = [0 20 280 360 380];
nz = length(zcoordinates);
a = zeros(n1,np);
rectarray=flip(rectarray);
max_row = size(rectarray,1);
max_column = size(rectarray,2);
  for ii = 1:1:np
    for irow = 1:1:max_row % bottom to top
       for icol = 1:1:max_column 
                im = rectarray(irow,icol);  
              if (location.x(ii)>=xcoordinates(1,3))&&(location.x(ii)<=xcoordinates(1,4))&&...
                 (location.y(ii)>=ycoordinates(3,1))&&(location.y(ii)<=ycoordinates(4,1))&&...
                 (location.z(ii)>=zcoordinates(1,3))&&(location.z(ii)<=zcoordinates(1,4))
                     a(:,ii) =  coefficient_a(xs(3));
                     break
              elseif (location.x(ii)>=xcoordinates(irow,icol))&&(location.x(ii)<=xcoordinates(irow,icol+1))&&...
                     (location.y(ii)>=ycoordinates(irow,icol))&&(location.y(ii)<=ycoordinates(irow+1,icol))&&...
                     (location.z(ii)>=zcoordinates(1,2))&&(location.z(ii)<=zcoordinates(1,4))
                     a(:,ii) =  coefficient_a(xs(im));
                      break    
              elseif(location.z(ii)>=zcoordinates(1,1))&&(location.z(ii)<=zcoordinates(1,2))
                     a(:,ii) =  coefficient_a(xs(4));
                     
              elseif(location.z(ii)>=zcoordinates(1,4))&&(location.z(ii)<=zcoordinates(1,5))&& im~=3
                     a(:,ii) =  coefficient_a(xs(4));
                     
              elseif(location.z(ii)>=zcoordinates(1,4))&&(location.z(ii)<=zcoordinates(1,5))&& im==3
                     a(:,ii) =  coefficient_a(xs(5));
              else
              end 
          
       end
    end
  end
end



















