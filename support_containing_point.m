function marked = support_containing_point(hspace,hmsh,points)

%INPUT:
%hspace: hierarchical space
%hmsh: hierarchcial mesh
%points: nx2 matrix containing coordinates of points that we want to use to select
%hierarchical basis functions
%OUTPUT:
%marked: cell array containing, for each level, the indices of the
%functions whose support contains at least one point in data(ind_data,:)

dx=hspace.space_of_level(1).degree(1);
dy=hspace.space_of_level(1).degree(2); 
ind=1;
mlev=1;
ndof_prev_levs=0;
clear marked
marked = cell (hspace.nlevels, 1); % cell-array che contiene, per ogni
% livello, gli indici delle funzioni i cui supporti sono da raffinare
for r=1:hspace.ndof
    while hspace.ndof_per_level(mlev)==0 && mlev<hspace.nlevels
        mlev=mlev+1;
    end
    rl=hspace.active{mlev}(r-ndof_prev_levs);
    [rx,ry]=ind2sub (hspace.space_of_level(mlev).ndof_dir, rl);  
    nx_m=hmsh.mesh_of_level(mlev).nel_dir(1);
    ny_m=hmsh.mesh_of_level(mlev).nel_dir(2);
    support=sp_get_cells(hspace.space_of_level(mlev),hmsh.mesh_of_level(mlev),(ry-1)*(nx_m+dx)+rx);
    [I,J]=ind2sub(hmsh.mesh_of_level(mlev).nel_dir, support); 
    %Indices of the knots which are the corners of the support (the indices you get are numbered starting from 1)
    I1_x=min(I);
    I2_x=max(I)+1;
    I1_y=min(J);
    I2_y=max(J)+1;    
    %Corresponding abscissae and ordinates
    x_mux=0+(I1_x-1)/nx_m;
    x_nux=0+(I2_x-1)/nx_m;
    y_muy=0+(I1_y-1)/ny_m;
    y_nuy=0+(I2_y-1)/ny_m;
    inx=find(points(:,1)>=x_mux & points(:,1)<=x_nux);
    iny=find(points(:,2)>=y_muy & points(:,2)<=y_nuy);    
    inner=intersect(inx,iny);

    if numel(inner)>0
        marked{mlev}(ind)=rl;
        ind=ind+1;
    end
    while r>=ndof_prev_levs+hspace.ndof_per_level(mlev) && mlev<hspace.nlevels
        ndof_prev_levs=ndof_prev_levs+hspace.ndof_per_level(mlev);
        mlev=mlev+1;
        ind=1;
    end
end
end

