%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to Compute Connectivity Matrix
% function C1 = ConnectivityMx(nf, np, faceMX)
%     C1 = sparse(nf, np);
%     for i = 1:nf
%         C1(i, faceMX(i,2)) = 1; 
%         C1(i, faceMX(i,3)) = -1;
%     end
% end
function C1 = ConnectivityMx(nf,np,faceMx)
    fi = (1:nf);
    p1 = faceMx(:,2)'; p2 = faceMx(:,3)'; % P1IDx   % P2 Idx
    C1=sparse([fi, fi], [p1, p2],[ones(1,nf), (-1)*ones(1,nf)],nf,np);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%