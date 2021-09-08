function OBJ1 = transformOBJCS6(varargin)





%OBJ1,figno,tissues,XYZ
if length(varargin)==3
    OBJ1 = varargin{1};
    theta = varargin{2};
    u = varargin{3};
    

    OBJ1.vertices = Quaternion3(theta,u,OBJ1.vertices);
    
     
    
else
     OBJ1 = varargin{1};
     flip = varargin{2};    
     
     
     if strcmp(flip,'flipX')==1
     OBJ1.vertices(:,1) = -OBJ1.vertices(:,1);
    
  
     elseif strcmp(flip,'flipY')==1
     OBJ1.vertices(:,2) = -OBJ1.vertices(:,2);         
         
                    
     elseif strcmp(flip,'flipZ')==1
     OBJ1.vertices(:,3) = -OBJ1.vertices(:,3);         
         
                    
     end
     
     
     
end