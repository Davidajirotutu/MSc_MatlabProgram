function power_array2 = power_rect_3D(result,rectarray,llcoordinate,...
                        lengtharray_x,lengtharray_y,xs,totalpower)
    
  % 1D Integration points Table provides int.point at interval [-1 1], so we Convert interval to[0,10] 
    Z(1,:) = -0.148874339*5*ones(1,7)+ 5.0; 
    Z(2,:) =  0.148874339*5*ones(1,7)+ 5.0; 
    Z(3,:) = -0.433395394*5*ones(1,7)+ 5.0;
    Z(4,:) =  0.433395394*5*ones(1,7)+ 5.0;
    Z(5,:) = -0.679409568*5*ones(1,7)+ 5.0;
    Z(6,:) =  0.679409568*5*ones(1,7)+ 5.0;
    Z(7,:) = -0.865063367*5*ones(1,7)+ 5.0;
    Z(8,:) =  0.865063367*5*ones(1,7)+ 5.0;
    Z(9,:) = -0.973906529*5*ones(1,7)+ 5.0;
    Z(10,:)=  0.973906529*5*ones(1,7)+ 5.0;
    % 1D Integration points Table provides weights at interval [-1 1], so we Convert interval to[0,10]
    C = 5.0*[0.295524225 0.295524225 0.269266719 0.269266719 0.219086363...
         0.219086363 0.149451349 0.149451349 0.066671344 0.066671344];
    ng = size(result.Eigenvectors,2);  % number of energy grooup : 
    np = size(result.Eigenvectors,1);  % number of result point :
    
    if mod(size(result.Eigenvectors,1),np) ~= 0  %                                      "size(u,1)/np" supposed to give us 2, if otherwise then there is an error
        error('u dimension mismaches the number of points')
    end 
    
    nx=length(lengtharray_x);
    ny=length(lengtharray_y);
    nz=length(Z);
    
    rectarray=flip(rectarray,1);
    power_array = zeros(ny-1,nx-1,nz);
    power_array1 = zeros(ny-1,nx-1);
    % 7 point integration
    a1 = 0.797426985353; a2= 0.059715871789;
    b1 = 0.101286507323; b2= 0.470142064105;
    w1 = 0.225; w2 = 0.125939180544; w3 = 0.132394152788;
    L=[1/3 1/3 1/3; ...
        a1 b1 b1; b1 a1 b1; b1 b1 a1; ...
        a2 b2 b2; b2 a2 b2; b2 b2 a2];
    W=[w1 w2 w2 w2 w3 w3 w3];
   % tmp2 = zeros(8,9,9,10);
  for iz = 1:1:nz  
      
    for iy = 1:1:ny-1 % Deduct 1 because of the algorithm that generate coordinate already added 1
        
       for ix=1:1:nx-1 % Deduct 1 because of the algorithm that generate coordinate already added 1
           im = rectarray(iy,ix);
           if im >0
               % rectangle element center coordinates
               cx =llcoordinate(1) + sum(lengtharray_x(1:ix)) - lengtharray_x(ix)/2;
               cy =llcoordinate(2) + sum(lengtharray_y(1:iy)) - lengtharray_y(iy)/2;
               % counterclockwise points on the rect edge. from x positive
               plist = [lengtharray_x(ix)/2 lengtharray_x(ix)/2 0 -lengtharray_x(ix)/2 ...
                       -lengtharray_x(ix)/2 -lengtharray_x(ix)/2 0 lengtharray_x(ix)/2; ...                % plist=[10 10 0 -10 -10 -10 0 10;
                        0 lengtharray_y(iy)/2 lengtharray_y(iy)/2 lengtharray_y(iy)/2 ...
                        0 -lengtharray_y(iy)/2 -lengtharray_y(iy)/2 -lengtharray_y(iy)/2];    %        0 10 10 10 0 -10 -10 -10]
               plist(1,:) = plist(1,:)+cx;
               plist(2,:) = plist(2,:)+cy;
               
               for ii = 1:1:8
                   if ii<8
                       trianglevertex = [cx plist(1,ii) plist(1,ii+1); ...  % vertex coorditates for each triangle in the reactor(i,i)
                                         cy plist(2,ii) plist(2,ii+1)];
                   else
                       trianglevertex = [cx plist(1,ii) plist(1,1); ...
                                         cy plist(2,ii) plist(2,1)]; %tri2grid
                   end

                   nintp=7;
                   IntPointCoordinates = zeros(2,nintp); % coordinates of all integration points ( 2 by 7)
                    for jj=1:1:nintp %calculate the integration point coordinates
                        IntPointCoordinates(1,jj)= dot(    trianglevertex(1,:), L(jj,:)   ); %x coordinate
                        IntPointCoordinates(2,jj)= dot(    trianglevertex(2,:), L(jj,:)   ); %y coordinate
                        
                    end
                    
                    IntPointCoordinates = [IntPointCoordinates(1,:);IntPointCoordinates(2,:);Z(iz,:)];
                    
                    tmp =zeros(7,ng);
                    for ig = 1:1:ng
                            tmp(:,ig) = interpolateSolution(result,IntPointCoordinates,ig,1);
                    end
                    
                    checknan = isnan(tmp);
                    
                    if sum(sum(checknan))<=0
                        
                         if sum(sum(tmp<0))>0
                           tmp = -tmp;
                         end
                       avgu=zeros(ng,1);
                      tmp = tmp';
                       for ig=1:1:ng
                            avgu(ig)=  dot(    tmp(ig,:),W(:)    ); 
                       end
                        tmp1 = dot(avgu(:),xs(im).vf) * polyarea(trianglevertex(1,:),trianglevertex(2,:) ) ; % integrated power in this triangle subdomain
                         power_array(iy,ix,iz) = power_array(iy,ix,iz) + tmp1;
                    end
               end
               
           end
            
       end  
       
    end
        power_array1 =( C(iz) * power_array(:,:,iz)) + power_array1;
  end
        power_array2 = (power_array1)/sum(sum(power_array1))*totalpower;
        power_array2 = flip(power_array2);
end