nc = num_clusters; %(number ofvertices on coarse grid); 
nv=num_vertices;
P0=sparse(nv,nc);

fo j = 1:num_clusters 
    P0(nodes_in_j_cluster,j)=ones(nodes_in_j_cluster,1);
end

Dinv=spdiags(1./diag(A),0,nv,nv);

tau=1.25;

B=(I-tau*Dinv*A);
P=B*P0;

%% two level
u=zeros(n,1);
r=b-A*u;
p=[1:n]';%'
e= GaussSeidel(A, r, p, 1e-10, 1);
u=u+e;
r=b-A*u;
%cg correction
Ac=P'*A*P;
ec=Ac\(P'*r);
u=u+P*ec;








