V = [F*K,-ones(size(F,1),sizes.nAlpha),zeros(size(F,1),sizes.nBeta);
    G*KW,zeros(size(G,1),sizes.nAlpha),-ones(size(G,1),sizes.nBeta);
    zeros(sizes.nAlpha+sizes.nBeta,sizes.nX),-eye(sizes.nAlpha+sizes.nBeta)];
v = ones(size(F,1)+size(G,1)+sizes.nAlpha+sizes.nBeta,1);

[V,v] = PreConditioning(V,v);

options = cplexoptimset('Display', 'off');

Vn =zeros(size(V));
vn = zeros(size(v));

Vpost = [];
vpost = [];

Vpre = V;
vpre = v;

iter = 1;
iterMax = 50;
while and(~altIsContained(Vpre,vpre,Vpost,vpost),iter<iterMax)
    
    if iter~=1
       [Vpre,vpre] = PreConditioning(Vpost,vpost);
        Vn = zeros(size(Vpre));
        vn = zeros(size(vpre));
    end
    
    for i = 1:size(Vpre,1)
    [x,fval,exitflag,output,lambda] = cplexlp(-Vpre(i,1:sizes.nX)*D,...
        G,ones(size(G,1),1),[],[],[],[],[],options);
    Vn(i,:) = [Vpre(i,1:sizes.nX)*(A+B*K),Vpre(i,sizes.nX+1),...
        Vpre(i,sizes.nX+sizes.nAlpha+1)+Vpre(i,1:sizes.nX)*x];
    vn(i) = vpre(i)-Vpre(i,1:sizes.nX)*D*x;
    end

[Vpost,vpost] = rowReduce([Vpre;Vn],[vpre;vn]);
fprintf('At iteration %2.0d the number of constraints was %3.0d.\n',iter,length(vpost))

iter = iter+1;
end

eqVec = [1;-1];

[Vpost,vpost] = PreConditioning(Vpost,vpost);
[Vout,vout] = rowReduce(Vpost,vpost);
[Vpost,vpost] = polyProject([Vout;[zeros(2),eqVec,zeros(2,1);zeros(2,3),eqVec]],...
    [vout;zeros(4,1)],2);

pInvSet = struct('A',Vpost,'b',vpost);

% Constant Terminal set computation. Terminal set computation is working.

V = [F*K;G*KW];
v = ones(size(V,1),1);

[V,v] = PreConditioning(V,v);

Vn =zeros(size(V));
vn = zeros(size(v));

Vpost = [];
vpost = [];

Vpre = V;
vpre = v;

iter = 1;
iterMax = 50;
while and(~altIsContained(Vpre,vpre,Vpost,vpost),iter<iterMax)
    
        if iter~=1
           [Vpre,vpre] = PreConditioning(Vpost,vpost);
            Vn = zeros(size(Vpre));
            vn = zeros(size(vpre));
        end

        for i = 1:size(Vpre,1)
        [x,fval,exitflag,output,lambda] = cplexlp(-Vpre(i,:)*D,G,ones(size(G,1),1));
        Vn(i,:) = Vpre(i,:)*(A+B*K);
        vn(i) = vpre(i)+fval;
        end
        
        [Vpost,vpost] = rowReduce([Vpre;Vn],[vpre;vn]);
        fprintf('At iteration %2.0d the number of constraints was %3.0d.\n',...
        iter,length(vpost))

    iter = iter+1;
end

[Vpost,vpost] = PreConditioning(Vpost,vpost);
cInvSet = struct('A',Vpost,'b',vpost);

figure(1)
plot(Polyhedron(cInvSet.A,cInvSet.b),'alpha',.2,'color',[0,1,0],...
    Polyhedron(pInvSet.A,pInvSet.b),'alpha',.2,'color',[1,0,0])
figure(gcf)