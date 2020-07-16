%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRE BLENDER, INITIAL SURFACE GENERATION AND SHOT PROJECTION CS5

%import data

path='C:\';
cd(path)
%import embryo data, as a matrix with x,y,z coordinates for each tissue (t)
rawdata.embryo=readtable('embryo_trans_E25A.csv');
%import shot data, as a matrix with x,y,z coordinates 
rawdata.shots=readtable('E25A_shots_trans_200319.csv');
%Optional add specific colours for tissues
colors = [[0 0.5 0];[0 0 1];[0.9290 0.6940 0.1250];[0.5 0.5 0.5];[1 0 0];[0.2,0.2,0.2];[1 0.5 0];[0 0.5 1]];
%tissue names, should match the order of tissues in your tissue collumn of
%embryo data
tissue_names=["Am","EmDisc","ExMes_stalk","ExMes_top","ExMes_bottom","SYS","Tb","VE"];
%matching blender order to erin original
%for centering manual imput of total size of sections
length_sys=(183-131)*12;

hold on
for tis = 0:(size(unique(rawdata.shots.t),1))-1
    
    %load tissue and shot coordinates for each tissue
      tissue.embryo=table2array(rawdata.embryo(rawdata.embryo.t==tis,5:7));    
     tissue.shots=table2array(rawdata.shots(rawdata.shots.t==tis,2:4));
  
    if isempty(tissue.embryo)==0
   
        if tissue_names(tis+1)=="ExMes_stalk"
            tissue_exmes_stalk.embryo=tissue.embryo;
                 tissue_exmes_stalk.surf=boundary(tissue.embryo(:,1),tissue.embryo(:,2),tissue.embryo(:,3),1);
        
        end
        
        if tissue_names(tis+1)=="ExMes_top"
            tissue_exmes_top.embryo=tissue.embryo;
                 tissue_exmes_top.surf=boundary(tissue.embryo(:,1),tissue.embryo(:,2),tissue.embryo(:,3),1);
            trisurf(tissue_exmes_top.surf,tissue_exmes_top.embryo(:,1),tissue_exmes_top.embryo(:,2),tissue_exmes_top.embryo(:,3),'FaceColor',colors(tis+1,:), 'edgecolor','none','FaceAlpha',0.5)
            csvwrite('ExMes_top_E25A_points.csv',tissue_exmes_top.embryo)
               csvwrite('ExMes_top_E25A_faces.csv',tissue_exmes_top.surf)
        end
        
        if tissue_names(tis+1)=="ExMes_bottom"
            tissue_exmes_bot.embryo=tissue.embryo;
                 tissue_exmes_bot.surf=boundary(tissue.embryo(:,1),tissue.embryo(:,2),tissue.embryo(:,3),1);
            trisurf(tissue_exmes_bot.surf,tissue_exmes_bot.embryo(:,1),tissue_exmes_bot.embryo(:,2),tissue_exmes_bot.embryo(:,3),'FaceColor',colors(tis+1,:), 'edgecolor','none','FaceAlpha',0.5)
            csvwrite('ExMes_bot_E25A_points.csv',tissue_exmes_bot.embryo)
               csvwrite('ExMes_bot_E25A_faces.csv',tissue_exmes_bot.surf)
        end
        
        if tissue_names(tis+1) == "Am" 

             tissue_am.embryo=tissue.embryo;
             shots_am=tissue.shots;
             length_amnion=(212-151)*12;
             
             amnion_max=max(tissue_am.embryo(:,3));
             amnion_min=min(tissue_am.embryo(:,3));
             amnion_end=amnion_max+(212-193)*12;
             proportion_covered=round(abs(amnion_max-amnion_min)/abs(amnion_end-amnion_min),2);
             
             
             %top and bottome centre of mass calc, to predict unsectioned
             %area
             amnion_max_tis=tissue_am.embryo(tissue_am.embryo(:,3)==amnion_max,:);
             minx = min(amnion_max_tis(:,1));
             maxx = max(amnion_max_tis(:,1));
             centx = (minx + maxx) / 2;
             miny = min(amnion_max_tis(:,2));
             maxy = max(amnion_max_tis(:,2));
             centy = (miny + maxy) / 2;
             center_mass_top=[centx,centy];
             
             
             amnion_min_tis=tissue_am.embryo(tissue_am.embryo(:,3)==amnion_min,:);
             minx = min(amnion_min_tis(:,1));
             maxx = max(amnion_min_tis(:,1));
             centx = (minx + maxx) / 2;
             miny = min(amnion_min_tis(:,2));
             maxy = max(amnion_min_tis(:,2));
             centy = (miny + maxy) / 2;
             center_mass_bot=[centx,centy];
             
             curve_seq=0:0.0001:10;
             prop_curve=(flip(sqrt(curve_seq)))/max(sqrt(curve_seq));
           
             tissue_am_embryo_plus=tissue_am.embryo;
             hold on
             for curve = 0:10-1
                 temp_prop=proportion_covered+(((1-proportion_covered)/10)*(curve+1));
                 
                 curve_subtrac=prop_curve(round(length(curve_seq)*temp_prop))/prop_curve(round(length(curve_seq)*proportion_covered));
            
                 new_am=amnion_max_tis;
                 vec=[center_mass_top,amnion_max]-new_am;
                 nvec = vec ./ sqrt(sum(vec .^ 2, 2)); 
                 
                                
                 dist=sqrt(sum(vec .^ 2, 2))*(1-curve_subtrac);
                 
                 final_points = new_am + dist.* nvec;
                 final_points(:,3)=(temp_prop*abs(amnion_end-amnion_min))+(amnion_min);
             
                 scatter3(final_points(:,1),final_points(:,2),final_points(:,3),...
                 40,'red','k','filled');
                 tissue_am_embryo_plus=vertcat(tissue_am_embryo_plus,final_points);
             end
            
             
             tissue_am.surf=boundary(tissue_am_embryo_plus(:,1),tissue_am_embryo_plus(:,2),tissue_am_embryo_plus(:,3),0.3);
            
        end
        if tissue_names(tis+1) == "EmDisc" 

             tissue_em.embryo=tissue.embryo;
             shots_em=tissue.shots;
             tissue_em.surf=boundary(tissue.embryo(:,1),tissue.embryo(:,2),tissue.embryo(:,3),0.1);
         
        end
    if tissue_names(tis+1) == "Tb" 
         tissue.embryo(:,1)=tissue.embryo(:,1)+100;
         tissue_Tb.embryo=tissue.embryo;
    
         tissue_Tb.surf=boundary(tissue.embryo(:,1),tissue.embryo(:,2),tissue.embryo(:,3),1);

  
         %Now see which Exmes in TB and project it up
            % find which which EMES in TB
            ineX_tb = inpolyhedron(tissue_Tb.surf,tissue_Tb.embryo,tissue_exmes_stalk.embryo);
            
             %project exmes onto trophoblast which overlapped
             [minProj,minDistances] = project2surf(tissue_Tb.surf,tissue_Tb.embryo,tissue_exmes_stalk.embryo(ineX_tb,:));
    
    
             [a,b]=min(minDistances,[],3);
 
            proj_points=zeros(3,size(tissue_exmes_stalk.embryo(ineX_tb,:),1));
            for pts = 1:size(tissue_exmes_stalk.embryo(ineX_tb,:),1)
                proj_points(:,pts,:)=minProj([1,2,3],pts,b(pts));
 
            end
            
                %then replace points
             tissue_exmes_stalk.embryo(ineX_tb,:)=proj_points.';
             
             
             %now same for amnion
              ineX_am = inpolyhedron(tissue_am.surf,tissue_am_embryo_plus,tissue_exmes_stalk.embryo);
                   
              tissue_exmes_stalk.embryo(ineX_am,:)=[];
              tissue_exmes_stalk.surf=boundary(tissue_exmes_stalk.embryo(:,1),tissue_exmes_stalk.embryo(:,2),tissue_exmes_stalk.embryo(:,3),0.8);

              csvwrite('Tb_E25A_points.csv',tissue_Tb.embryo)
               csvwrite('Tb_E25A_faces.csv',tissue_Tb.surf)
               
               csvwrite('ExMes_stalk_E25A_points.csv',tissue_exmes_stalk.embryo)
               csvwrite('ExMes_stalk_E25A_faces.csv',tissue_exmes_stalk.surf)
                trisurf(tissue_Tb.surf,tissue_Tb.embryo(:,1),tissue_Tb.embryo(:,2),tissue_Tb.embryo(:,3),'FaceColor',colors(tis+1,:), 'edgecolor','none','FaceAlpha',0.5)             
      
    end
    if tissue_names(tis+1) == "SYS" 
             tissue_sys.embryo=tissue.embryo;
             shots_sys=tissue.shots;
             length_sys=(183-131)*12;
             
             sys_max=max(tissue_sys.embryo(:,3));
             sys_min=min(tissue_sys.embryo(:,3));
             sys_end=sys_min-(151-131)*12;
             proportion_covered=round(abs(sys_max-sys_min)/abs(sys_max-sys_end),2);
             
             
   %top and bottome centre of mass calc, to predict unsectioned
             %area
             sys_max_tis=tissue_sys.embryo(tissue_sys.embryo(:,3)==sys_max,:);
             minx = min(sys_max_tis(:,1));
             maxx = max(sys_max_tis(:,1));
             centx = (minx + maxx) / 2;
             miny = min(sys_max_tis(:,2));
             maxy = max(sys_max_tis(:,2));
             centy = (miny + maxy) / 2;
             center_mass_top=[centx,centy];
             
             
             sys_min_tis=tissue_sys.embryo(tissue_sys.embryo(:,3)==sys_min,:);
             minx = min(sys_min_tis(:,1));
             maxx = max(sys_min_tis(:,1));
             centx = (minx + maxx) / 2;
             miny = min(sys_min_tis(:,2));
             maxy = max(sys_min_tis(:,2));
             centy = (miny + maxy) / 2;
             center_mass_bot=[centx,centy];
             
             curve_seq=0:0.0001:10;
             prop_curve=(flip(sqrt(curve_seq)))/max(sqrt(curve_seq));
           
                tissue_sys_embryo_plus=tissue_sys.embryo;
             hold on
             for curve = 0:10-1
                 temp_prop=proportion_covered+(((1-proportion_covered)/10)*(curve+1));
                 
                 curve_subtrac=prop_curve(round(length(curve_seq)*temp_prop))/prop_curve(round(length(curve_seq)*proportion_covered));
            
                 new_sys=sys_min_tis;
                 vec=[center_mass_bot,sys_min]-new_sys;
                 nvec = vec ./ sqrt(sum(vec .^ 2, 2)); 
                 
                                
                 dist=sqrt(sum(vec .^ 2, 2))*(1-curve_subtrac);
                 
                 final_points = new_sys + dist.* nvec;
                 final_points(:,3)=(temp_prop*(sys_end-sys_max))+(sys_max);
              final_points(:,1)=final_points(:,1)+250;
              final_points(:,2)=final_points(:,2)-250;
                 scatter3(final_points(:,1),final_points(:,2),final_points(:,3),...
                 40,'red','k','filled');
                 tissue_sys_embryo_plus=vertcat(tissue_sys_embryo_plus,final_points);
             end
            
             
             tissue_sys.surf=boundary(tissue_sys_embryo_plus(:,1),tissue_sys_embryo_plus(:,2),tissue_sys_embryo_plus(:,3),0.1);
      
             trisurf(tissue_sys.surf,tissue_sys_embryo_plus(:,1),tissue_sys_embryo_plus(:,2),tissue_sys_embryo_plus(:,3),'FaceColor',colors(tis+1,:), 'edgecolor','none','FaceAlpha',0.5)
     
             %check amnion points in yolk sac and move out
            in_sys = inpolyhedron(tissue_sys.surf,tissue_sys_embryo_plus,tissue_am_embryo_plus);
              %then we project those points onto the surface of the
              %trophoblast
           [minProj,minDistances] = project2surf(tissue_sys.surf, tissue_sys_embryo_plus,tissue_am_embryo_plus(in_sys,:));
           [a,b]=min(minDistances,[],3);

           proj_points=zeros(3,size(tissue_am_embryo_plus(in_sys,:),1));
               for pts = 1:size(tissue_am_embryo_plus(in_sys,:),1)
                   proj_points(:,pts,:)=minProj([1,2,3],pts,b(pts));

               end
               new_amSYS=proj_points;
               %then replace points
               tissue_am_embryo_plus(in_sys,:)=new_amSYS.';
               tissue_am.surf=boundary(tissue_am_embryo_plus(:,1),tissue_am_embryo_plus(:,2),tissue_am_embryo_plus(:,3),0.3);
      
               
         %Now do the reverse
             in_am = inpolyhedron(tissue_am.surf,tissue_am_embryo_plus,tissue_sys_embryo_plus);
              %then we project those points onto the surface of the
              %trophoblast
           [minProj,minDistances] = project2surf(tissue_am.surf, tissue_am_embryo_plus,tissue_sys_embryo_plus(in_am,:));
           [a,b]=min(minDistances,[],3);

           proj_points=zeros(3,size(tissue_sys_embryo_plus(in_am,:),1));
               for pts = 1:size(tissue_sys_embryo_plus(in_am,:),1)
                   proj_points(:,pts,:)=minProj([1,2,3],pts,b(pts));

               end
               new_sysAM=proj_points;
               %then replace points
               tissue_sys_embryo_plus(in_am,:)=[];
            tissue_sys.surf=boundary(tissue_sys_embryo_plus(:,1),tissue_sys_embryo_plus(:,2),tissue_sys_embryo_plus(:,3),0.4);
         
      %Now check for points in the epiblast in both tissues
         
                    in_epi_sys = inpolyhedron(tissue_em.surf,tissue_em.embryo,tissue_sys_embryo_plus);
              %then we project those points onto the surface of the
              %SYS
           [minProj,minDistances] = project2surf(tissue_em.surf, tissue_em.embryo,tissue_sys_embryo_plus(in_epi_sys,:));
           [a,b]=min(minDistances,[],3);

           proj_points=zeros(3,size(tissue_sys_embryo_plus(in_epi_sys,:),1));
               for pts = 1:size(tissue_sys_embryo_plus(in_epi_sys,:),1)
                   proj_points(:,pts,:)=minProj([1,2,3],pts,b(pts));

               end
               new_epiSYS=proj_points;
               %then replace points
               tissue_sys_embryo_plus(in_epi_sys,:)=new_epiSYS.';
               tissue_sys.surf=boundary(tissue_sys_embryo_plus(:,1),tissue_sys_embryo_plus(:,2),tissue_sys_embryo_plus(:,3),0.8);
         
      
         in_epi_am = inpolyhedron(tissue_em.surf,tissue_em.embryo,tissue_am_embryo_plus);
     %then we check amnion
           [minProj,minDistances] = project2surf(tissue_em.surf, tissue_em.embryo,tissue_am_embryo_plus(in_epi_am,:));
           [a,b]=min(minDistances,[],3);

           proj_points=zeros(3,size(tissue_am_embryo_plus(in_epi_am,:),1));
               for pts = 1:size(tissue_am_embryo_plus(in_epi_am,:),1)
                   proj_points(:,pts,:)=minProj([1,2,3],pts,b(pts));

               end
               new_epiAM=proj_points;
               %then replace points
               tissue_am_embryo_plus(in_epi_am,:)=new_epiAM.';
               
    
                          %regenerate surfaces
            tissue_am.surf=boundary(tissue_am_embryo_plus(:,1),tissue_am_embryo_plus(:,2),tissue_am_embryo_plus(:,3),0.6);
            trisurf(tissue_am.surf,tissue_am_embryo_plus(:,1),tissue_am_embryo_plus(:,2),tissue_am_embryo_plus(:,3),'FaceColor',colors(tis+1,:), 'edgecolor','none','FaceAlpha',0.5)
               
            %now do the same for the epiblast 
            in_am_epi = inpolyhedron(tissue_am.surf,tissue_am_embryo_plus,tissue_em.embryo);
    
                      tissue_em.embryo(in_am_epi,:)=[];
               
    
                          %regenerate surfaces
            tissue_em.surf=boundary(tissue_em.embryo(:,1),tissue_em.embryo(:,2),tissue_em.embryo(:,3),1);
            %trisurf(tissue_em.surf,tissue_em.embryo(:,1),tissue_em.embryo(:,2),tissue_em.embryo(:,3),'FaceColor',colors(2,:), 'edgecolor','none','FaceAlpha',0.5)
               
            in_sys_epi = inpolyhedron(tissue_sys.surf,tissue_sys_embryo_plus,tissue_em.embryo);
     
               tissue_sys_embryo_plus(in_sys_epi,:)=[];
               
                              %regenerate surfaces
             tissue_em.surf=boundary(tissue_em.embryo(:,1),tissue_em.embryo(:,2),tissue_em.embryo(:,3),1);
            trisurf(tissue_em.surf,tissue_em.embryo(:,1),tissue_em.embryo(:,2),tissue_em.embryo(:,3),'FaceColor',colors(tis+1,:), 'edgecolor','none','FaceAlpha',0.7)
                   
            csvwrite('EmDisc_E25A_points.csv',tissue_em.embryo)
               csvwrite('EmDisc_E25A_faces.csv',tissue_em.surf)
         
            csvwrite('SYS_E25A_points.csv',tissue_sys_embryo_plus)
               csvwrite('SYS_E25A_faces.csv',tissue_sys.surf)
               
            csvwrite('Am_E25A_points.csv',tissue_am_embryo_plus)
               csvwrite('Am_E25A_faces.csv',tissue_am.surf)
                
    end
    
   
    if tissue_names(tis+1)=="VE"
         
         tissue_VE.embryo=tissue.embryo;
        
                             
               [minProj,minDistances] = project2surf(tissue_sys.surf, tissue_sys_embryo_plus,tissue_VE.embryo);

               [a,b]=min(minDistances,[],3);

               proj_points=zeros(3,size(tissue_VE.embryo,1));
               for pts = 1:size(tissue_VE.embryo,1)
                   proj_points(:,pts,:)=minProj([1,2,3],pts,b(pts));

               end
               
              tissue_VE_embryo_plus=proj_points;
               
         tissue_VE.surf=boundary(tissue_VE_embryo_plus(:,1),tissue_VE_embryo_plus(:,2),tissue_VE_embryo_plus(:,3),1);
        
         shots_VE=tissue.shots;
         csvwrite('VE_E25A_points.csv',tissue_VE_embryo_plus)
               csvwrite('VE_E25A_faces.csv',tissue_VE.surf)
     end
   

   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POST BLENDER

% Am_ 0
% EmDisc_ 1
% ExMes_stalk_ 2
% ExMes_top _3
% ExMes_bottom _ 4
% SYS_ 5
% Tb_ 6
% VE_ 7
%load blender file, require read_wobj_ds package
blend_file=read_wobj_ds;
tis_nums=[1,4,8,12,16,20];
smoothedDT=blend_file.vertices;
colors = [[0 0.5 0];[0 0 1];[0.9290 0.6940 0.1250];[0.5 0.5 0.5];[1 0 0];[0.2,0.2,0.2];[1 0.5 0];[0 0.5 1]];
%check these always as each obj file can be diff, but hard code rearrage
%tissue names
tissue_names=["Tb","SYS","VE","Am","EmDisc","ExMes_stalk","ExMes"];
tissue_names_old=["Am","EmDisc","ExMes_stalk","ExMes","ExMes","SYS","Tb","VE"];
%rearrange shots to match order
shot_tissues=strings(size(rawdata.shots.t,1),1);
for tis = 1:(size(rawdata.shots.t,1))
    shot_tissues(tis,1)=tissue_names_old(rawdata.shots.t(tis)+1);
end
prism_z=[180,204,228,240,252,300];
min_troph=min(blend_file.objects(4).data.vertices(:,3));

%loop through tissues
points_table_final=cell(0,4);
for tis = 0:(size(tissue_names,2))-1
         %load shots and tissue coordinates for tissue being interrogated
    if tissue_names(tis+1) == "VE" 
        tissue.shots=table2array(rawdata.shots(shot_tissues=="VE",2:4));
        tissue.names=rawdata.shots(shot_tissues=="VE",1);
         hull = blend_file.objects(4*(tis+1)).data.vertices;
    elseif tissue_names(tis+1) == "Tb" 
        tissue.shots=table2array(rawdata.shots(shot_tissues=="Tb",2:4));
        tissue.names=rawdata.shots(shot_tissues=="Tb",1);
         hull = blend_file.objects(4*(tis+1)).data.vertices;
    elseif tissue_names(tis+1) == "ExMes" 
         tissue.shots=table2array(rawdata.shots(shot_tissues=="ExMes" ,2:4));
        tissue.names=rawdata.shots(shot_tissues=="ExMes",1);
         hull = blend_file.objects(34).data.vertices;
    elseif tissue_names(tis+1) == "ExMes_stalk" 
         tissue.shots=table2array(rawdata.shots(shot_tissues=="ExMes_stalk" ,2:4)); 
        tissue.names=rawdata.shots(shot_tissues=="ExMes_stalk",1);
     hull = blend_file.objects(30).data.vertices;
    elseif tissue_names(tis+1) == "Am" 
         tissue.shots=table2array(rawdata.shots(shot_tissues=="Am",2:4));
        tissue.names=rawdata.shots(shot_tissues=="Am",1);
         hull = blend_file.objects(4*(tis+1)).data.vertices;
    elseif tissue_names(tis+1) == "SYS" 
        tissue.shots=table2array(rawdata.shots(shot_tissues=="SYS",2:4));
        tissue.names=rawdata.shots(shot_tissues=="SYS",1);
         hull = blend_file.objects(4*(tis+1)).data.vertices;
    elseif tissue_names(tis+1) == "EmDisc" 
        tissue.shots=table2array(rawdata.shots(shot_tissues=="EmDisc",2:4));
        epi_verts= smoothedDT(unique(blend_file.objects(4*(tis+1)).data.vertices),:);
        z_epi=unique(tissue.shots(:,3));
        %manual proportioning base on measurements of IF, subtract to make first point match
        %proportion, and then multiply to get last point to match, this is
        %
        
        prop_161=[single(0.1);single(0.95)];
        prop_166=[single(0.5);single(0.7)];
        prop_169=[single(0.1);single(0.80)];
        prop_174=[single(0.1);single(0.80)];
        prop_180=[single(0.02);single(0.85)];
        prop_184=[single(0.58);single(0.62)];
        prop=horzcat(prop_161,prop_166,prop_169,prop_174,prop_180,prop_184);
        %look at 169
        %184 make cluster tighter
        for epis=1:size(z_epi,1)
            temp_shots=tissue.shots(tissue.shots(:,3)==z_epi(epis),:);
            temp_tis=epi_verts((z_epi(epis)+20)>=epi_verts(:,3)& epi_verts(:,3)>=(z_epi(epis)-20),:);
            max_x=max(temp_tis(:,1));
            min_x=min(temp_tis(:,1));
            shift_x=max(temp_shots(:,1))-((prop(1,epis)*(min_x-max_x))+max_x);
            multiplier_x=((prop(2,epis)*(min_x-max_x))+max_x)/(min(temp_shots(:,1)-shift_x));
            new_shots=temp_shots;
            new_shots(:,1)=(new_shots(:,1)-shift_x)*multiplier_x;
            tissue.shots(tissue.shots(:,3)==z_epi(epis),1)=new_shots(:,1);
        end
       
        tissue.names=rawdata.shots(shot_tissues=="EmDisc",1);
         hull = blend_file.objects(4*(tis+1)).data.vertices;
    end

    hold on
%plot surfaces
    trisurf(hull,smoothedDT(:,1),smoothedDT(:,2),smoothedDT(:,3),...
    'facealpha',.3,'LineWidth', 0.01,'Facecolor',colors(tis+1,:), 'edgecolor','none','FaceAlpha',0.5)
%project points onto new surfaces
   [minProj,minDistances] = project2surf(hull, smoothedDT,tissue.shots);
   [a,b]=min(minDistances,[],3);

   proj_points=zeros(3,size(tissue.shots,1));
   for pts = 1:size(tissue.shots,1)
       proj_points(:,pts,:)=minProj([1,2,3],pts,b(pts));
       
   end
    tisse_name_export    = cell(size(permute(proj_points,[2,1]),1), 1);
    tisse_name_export(:) = {tissue_names(tis+1)};
    
    points_table=horzcat(num2cell(permute(proj_points,[2,1])),table2array(tissue.names),tisse_name_export);
    points_table_final=vertcat(points_table_final,points_table);
     scatter3(proj_points(1,:),proj_points(2,:),proj_points(3,:),...
                    40,'k','filled');
      csvwrite(strcat(tissue_names(tis+1),'_BLEND_points_high_final.csv'),smoothedDT)
      csvwrite(strcat(tissue_names(tis+1),'_BLEND_faces_high_final.csv'),hull)
     
end
final_shots_table=cell2table(points_table_final);
final_shots_table.Properties.VariableNames = ["x","y","z","shot","tissue"];
writetable(final_shots_table,'Proj_shots_BLEND_high.csv')

% Auxillary Functions======================================================    
% =========================================================================    

function [projX0,minDistances] = project2surf(vertices,points,shot)
%vertices are the points on the convex hull
%points are the points in the triangulation
%shots are the shots of interest

    % Get the number of faces and the number cells
    nPoints = size(points,1);
    nFaces = size(vertices,1);
    nShots = size(shot,1);
    
    % Build array of faces
    faces = cell(1,nFaces);    
    for i = 1:nFaces
        faces{i} = points(vertices(i,:),:);
    end
    
    % Calcualte the distace of each cell to each face and store the data in
    % distances -----------------------------------------------------------
    %
    % The distance of a point (x0) from a finite plane S, given three
    % points (x1,x2,x3) in the plane is calculated by:   
    % (1) Calculate the normal vector to the plane,
    %    n = cross(x2-x1,x3-x1)/norm(cross(x2-x1,x3-x1)).
    % (2) Project the difference of x0 and one of the point onto the normal
    % vector,
    %    distance = dot(n,x0-x1).
    % (3) Determine if the proj(x0) onto the coplanar surface of S is
    % within the domain of S
    % (4) If not, Find the closet point of the boundary of S.

    % Calculate the normal vector to each faces
    normVecs = cellfun(@(f) cross(f(1,:)-f(2,:),f(1,:)-f(3,:))/...
        norm( cross(f(1,:)-f(2,:),f(1,:)-f(3,:))),faces,...
        'UniformOutput',false);
   
        % Project the difference vector onto the normal vector
        X0 = repmat(shot',1,1,nFaces);
        
        X1 = repmat(points(vertices(:,1),:)',1,1,nShots);
        X1 = permute(X1, [1 3 2]);
        diffMat = X0-X1;

        normMat =  repmat(cat(1,normVecs{:})',1,1,nShots);
        normMat = permute(normMat, [1 3 2]);
        minSurfDistances = dot(normMat,diffMat);


        % Determine if the proj(X0) is in the triangluar domain----------------

        % calculate proj(X0) 
        %projected point onto the infinite plane
        surfProjX0 = X0 - normMat.*repmat(minSurfDistances,3,1,1);
        minSurfDistances = abs(dot(normMat,diffMat));

        X2 = repmat(points(vertices(:,2),:)',1,1,nShots);
        X2 = permute(X2, [1 3 2]);

        X3 = repmat(points(vertices(:,3),:)',1,1,nShots);
        X3 = permute(X3, [1 3 2]);

        u = X2 - X1;
        v = X3 - X1;
        w = X0 - X1;
        %normal 
        n = cross(u,v);

        gamma = dot(cross(u,w),n)./dot(n,n); %parallel test
        beta  = dot(cross(w,v),n)./dot(n,n);
        alpha = 1 - gamma - beta;

        T =   permute(gamma,[ 2 3 1])<=1 & permute(gamma,[ 2 3 1])>=0 ...
            & permute(beta, [ 2 3 1])<=1 & permute(beta, [ 2 3 1])>=0 ...
            & permute(alpha,[ 2 3 1])<=1 & permute(alpha,[ 2 3 1])>=0;
        
        surfProjX0(:,~T) = 0;
        minSurfDistances(:,~T) = 0;

        % Calculate the projection on the point to the line--------------------
        [p12, d12] = minDistPointOnSegment(X1,X2,X0);
        [p13, d13] = minDistPointOnSegment(X1,X3,X0);     
        [p23, d23] = minDistPointOnSegment(X2,X3,X0);


        D(:,:,1) = permute(d12,[2 3 1]);
        D(:,:,2) = permute(d13,[2 3 1]);
        D(:,:,3) = permute(d23,[2 3 1]);

        [minEdgeDistance,minEdgeIndex] = min(D,[],3);

        minEdgeDistance = permute(minEdgeDistance, [3 1 2] );

        minEdgeIndex = permute(minEdgeIndex, [3 1 2] );  

        p12(:,minEdgeIndex~=1) = 0;
        p13(:,minEdgeIndex~=2) = 0;
        p23(:,minEdgeIndex~=3) = 0;

        edgeProjX0 = p12+p13+p23;

        edgeProjX0(:,T) = 0;
        minEdgeDistance(:,T) = 0;

        projX0 = surfProjX0 + edgeProjX0;
        minDistances = minSurfDistances + minEdgeDistance;
    
end

function [projX0,minDistances] = distance2surf(vertices, points)

    % Get the number of faces and the number cells
    nPoints = size(points,1);
    nFaces = size(vertices,1);
    
    % Build array of faces
    faces = cell(1,nFaces);    
    for i = 1:nFaces
        faces{i} = points(vertices(i,:),:);
    end
    
    % Calcualte the distace of each cell to each face and store the data in
    % distances -----------------------------------------------------------
    %
    % The distance of a point (x0) from a finite plane S, given three
    % points (x1,x2,x3) in the plane is calculated by:   
    % (1) Calculate the normal vector to the plane,
    %    n = cross(x2-x1,x3-x1)/norm(cross(x2-x1,x3-x1)).
    % (2) Project the difference of x0 and one of the point onto the normal
    % vector,
    %    distance = dot(n,x0-x1).
    % (3) Determine if the proj(x0) onto the coplanar surface of S is
    % within the domain of S
    % (4) If not, Find the closet point of the bounary of S.

    % Calculate the normal vector to each faces
    normVecs = cellfun(@(f) cross(f(2,:)-f(1,:),f(3,:)-f(1,:))/...
        norm( cross(f(2,:)-f(1,:),f(3,:)-f(1,:))),faces,...
        'UniformOutput',false);
    
    % Project the difference vector onto the normal vector
    X0 = repmat(points',1,1,nFaces);
    X1 = repmat(points(vertices(:,1),:)',1,1,nPoints);
    X1 = permute(X1, [1 3 2]);
    diffMat = X0-X1;

    normMat =  repmat(cat(1,normVecs{:})',1,1,nPoints);
    normMat = permute(normMat, [1 3 2]);
    minSurfDistances = dot(normMat,diffMat);
    
    
    % Determine if the proj(X0) is in the triangluar domain----------------
    
    % calculate proj(X0) onto finite plane S
    surfProjX0 = X0 - normMat.*repmat(minSurfDistances,3,1,1);
    
    
    X2 = repmat(points(vertices(:,2),:)',1,1,nPoints);
    X2 = permute(X2, [1 3 2]);

    X3 = repmat(points(vertices(:,3),:)',1,1,nPoints);
    X3 = permute(X3, [1 3 2]);

    u = X2 - X1;
    v = X3 - X1;
    w = X0 - X1;
    
    n = cross(u,v);
    
    gamma = dot(cross(u,w),n)./dot(n,n);
    beta  = dot(cross(w,v),n)./dot(n,n);
    alpha = 1 - gamma - beta;
    
    T =   permute(gamma,[ 2 3 1])<=1 & permute(gamma,[ 2 3 1])>=0 ...
        & permute(beta, [ 2 3 1])<=1 & permute(beta, [ 2 3 1])>=0 ...
        & permute(alpha,[ 2 3 1])<=1 & permute(alpha,[ 2 3 1])>=0;
    
    surfProjX0(:,~T) = 0;
    minSurfDistances(:,~T) = 0;
    
    % Calculate the projection on the point to the line--------------------
    [p12, d12] = minDistPointOnSegment(X1,X2,X0);
    [p13, d13] = minDistPointOnSegment(X1,X3,X0);     
    [p23, d23] = minDistPointOnSegment(X2,X3,X0);
      
        
    D(:,:,1) = permute(d12,[2 3 1]);
    D(:,:,2) = permute(d13,[2 3 1]);
    D(:,:,3) = permute(d23,[2 3 1]);
    
    [minEdgeDistance,minEdgeIndex] = min(D,[],3);
    
    minEdgeDistance = permute(minEdgeDistance, [3 1 2] );
    
    minEdgeIndex = permute(minEdgeIndex, [3 1 2] );  
    
    p12(:,minEdgeIndex~=1) = 0;
    p13(:,minEdgeIndex~=2) = 0;
    p23(:,minEdgeIndex~=3) = 0;
    
    edgeProjX0 = p12+p13+p23;
    
    edgeProjX0(:,T) = 0;
    minEdgeDistance(:,T) = 0;
    
    projX0 = surfProjX0 + edgeProjX0;
    minDistances = minSurfDistances + minEdgeDistance;

end


function [point, distance] = minDistPointOnSegment(a,b,p)
    
    S = sqrt(sum((a-b).^2,1));
    A = sqrt(sum((a-p).^2,1));
    B = sqrt(sum((b-p).^2,1));
    
    % Calculate projection distance and vector
    ba = b - a; 
    pa = p - a;
    pb = p - b;
    
    normAB = cross(a,b)./ repmat(sqrt(sum(cross(a,b).^2,1)),3,1,1);
    dProjP = sqrt(sum(cross(ba,pa).^2,1))./ sqrt(sum(ba.^2,1));
    projP = a + repmat(dot(pa,ba)./dot(ba,ba),3,1,1).*ba;
    
    % Ensure projected point is on the line 
    orthoganalityTest = dot(normAB,projP-a)<1e-15;
    assert(all(orthoganalityTest(:)), 'projection not normal')
    
    % Determing if the projected point is within the line segment
    onOrOff = A.^2 + B.^2 - 2.*(dProjP.^2) < S.^2;
    
    projP(:,~onOrOff) = 0;
    dProjP(:,~onOrOff) = 0;
    % get the minEndPoint and Distance
    endPointDistance(:,:,1) = permute(sqrt(sum(pa.^2,1)),[2 3 1]);
    endPointDistance(:,:,2) = permute(sqrt(sum(pb.^2,1)),[2 3 1]);
    
    [minEndPointDistance,minEndPointIndex] =  min(endPointDistance,[],3);
    
    minEndPointDistance = permute(minEndPointDistance, [3 1 2] );
    minEndPointIndex = permute(minEndPointIndex, [3 1 2] );
    
    a(repmat(minEndPointIndex==2,3,1,1))=0; 
    b(repmat(minEndPointIndex==1,3,1,1))=0;
    
    minEndPoint = a+b;
    minEndPoint(:,onOrOff) = 0;
    
    minEndPointDistance(:,onOrOff) = 0;
    
    % Now us onOrOff to get the correct minimem distance to the line
    % segment
    point = projP +  minEndPoint;
    distance = dProjP+ minEndPointDistance;

end
