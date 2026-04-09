%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = makeDiffusiveConductance(faceMx,ptCoordMx)
    nf = size(faceMx,1);
    faces_vec = ptCoordMx(faceMx(:,3),:) - ptCoordMx(faceMx(:,2),:);
    faces_x = find(faces_vec(:,1)>0); hx = vecnorm(faces_vec(faces_x(1),:),2,2);
    faces_y = find(faces_vec(:,2)>0); hy = vecnorm(faces_vec(faces_y(1),:),2,2);
    faces_z = find(faces_vec(:,3)>0); hz = vecnorm(faces_vec(faces_z(1),:),2,2);
    conductance = nan(nf,1);
    conductance(faces_x) = hy*hz/hx;
    conductance(faces_y) = hx*hz/hy;
    conductance(faces_z) = hx*hy/hz;
    W = spdiags(conductance,0,nf,nf);
end
