function OBJ1 = transformOBJ(varargin)

%OBJ1,figno,tissues,XYZ
if length(varargin)==3
    OBJ1 = varargin{1};
    theta = varargin{2};
    u = varargin{3};
    

    try
    OBJ1.vertices = Quaternion3(theta,u,OBJ1.vertices);
    catch
        keyboard
        OBJ1.vertices = Quaternion3(theta,u,OBJ1.cleanX);
    end
         
    
else
     OBJ1 = varargin{1};
     flip = varargin{2};    
     
     
     if strcmp(flip,'flipX')==1
         try
     OBJ1.vertices(:,1) = -OBJ1.vertices(:,1);
         catch
                  OBJ1.cleanX(:,1) = -OBJ1.cleanX(:,1);

         end
    
  
     elseif strcmp(flip,'flipY')==1
         
         try
     OBJ1.vertices(:,2) = -OBJ1.vertices(:,2);
         catch
                  OBJ1.cleanX(:,2) = -OBJ1.cleanX(:,2);

         end
         
                    
     elseif strcmp(flip,'flipZ')==1
         try
     OBJ1.vertices(:,3) = -OBJ1.vertices(:,3);      
         catch
                  OBJ1.cleanX(:,3) = -OBJ1.cleanX(:,3);      

         end
         
                    
     end
     
     
     
end