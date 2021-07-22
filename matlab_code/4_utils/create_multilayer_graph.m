function [G3] = create_multilayer_graph(G,T,param)

% current: unweighted, directed or undirected, periodic or nonperiodic

N=G.N;
G2=G;

if param.weighted==0
    G2.W=G2.W>0;
    G2.d=sum(G2.W);
    G2.d=G2.d(:);
    G2.L=diag(G2.d)-G2.W;
end

G3=struct;
G3.N=G2.N*T;
G3.Ne=G2.Ne*T+N*T;
G3.directed=param.directed;

if param.directed
    G3.d=repmat(G2.d,T,1)+2; % check: +1 or +2 for directed case? weights of new edges? 
    % the function gsp_jtv_graph creates unweighted edges
    pos=diag(ones(T-1,1),1);
    if param.periodic==1
        pos(T,1)=1;
    end
else
    G3.d=repmat(G2.d,T,1)+2;
    pos=diag(ones(T-1,1),1)+diag(ones(T-1,1),-1);
    if param.periodic==1
        pos(1,T)=1;
        pos(T,1)=1;
    end
end

% give an option for nonperiodic -> path graph instead of ring graph in
% time

G3.W=kron(pos,eye(N)); % change eye(N) to add weights for edges between layer
G3.W=G3.W-diag(diag(G3.W))+kron(eye(T),G2.W);
G3.L=diag(G3.d)-G3.W;
G3.coords=[repmat(G2.coords,T,1) vec(repmat(1:T,N,1))/T]; % divide by T at 
% the end to get z coordinates between 0 and 1, similar scale with x,y
horizontal_shift=0; % make layers offset to see the ring shape between them
G3.coords(N+1:2*N,1)=G3.coords(N+1:2*N,1)+horizontal_shift*ones(N,1); 
% TODO: adapt this offset for layer T>3 case
G3.plotting.vertex_size=G2.plotting.vertex_size;

end
