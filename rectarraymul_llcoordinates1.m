%  This function generate the Rectangle Array Lower Left Cordinate 
function [xcoordinates,ycoordinates] = rectarraymul_llcoordinates(rectarray,lengtharray_x,lengtharray_y,llcoordinate)
    max_column = size(lengtharray_x,2);
    max_row = size(lengtharray_y,2);
    
    for irow = 1:1:max_row % bottom to top
        for icol = 1:1:max_column % left to right
           if irow==1 && icol==1
               rect_x = llcoordinate(1,1) + (icol-1)*(lengtharray_x(icol)-lengtharray_x(icol))/2;
               rect_y = llcoordinate(1,2) + (irow-1)*lengtharray_y(irow);   
           elseif irow==1
               rect_x = llcoordinate(1,1) + (icol-1)*lengtharray_x(icol)-lengtharray_x(icol)/2;
               rect_y = llcoordinate(1,2) + (irow-1)*lengtharray_y(irow);  
           elseif icol==1
           rect_x = llcoordinate(1,1) + (icol-1)*(lengtharray_x(icol)-lengtharray_x(icol)/2);
           rect_y = llcoordinate(1,2) + (irow-1)*lengtharray_y(irow)-(lengtharray_y(irow)/2);    
                    
           else 
           rect_x = llcoordinate(1,1) + (icol-1)*lengtharray_x(icol)-(lengtharray_x(icol)/2);
           rect_y = llcoordinate(1,2) + (irow-1)*lengtharray_y(irow)-(lengtharray_y(irow)/2);  
           end
           domaincoordinates(irow,icol,1) = rect_x;
           domaincoordinates(irow,icol,2) = rect_y;
                 
        end
    end
xcoordinates = domaincoordinates(:,:,1);
ycoordinates = domaincoordinates(:,:,2);
rectcoordinates= domaincoordinates;
end

 
           