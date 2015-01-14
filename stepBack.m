function [Aout,bout] = stepBack(Ain,bin,s)

Vn = zeros(size(Ain,1),s.nX+s.nAlpha+s.nBeta+s.nU);
vn = zeros(size(bin));

for i = 1:size(Ain,1)
    [~,fval,~,~,~] = cplexlp(-Ain(i,1:s.nX)*s.D,s.G,ones(size(s.G,1),1));
    Vn(i,:) = [Ain(i,1:s.nX)*s.A,Ain(i,s.nX+1),Ain(i,s.nX+s.nAlpha+1)-fval,Ain(i,1:s.nX)*s.B];
    vn(i) = bin(i)+fval;
end

V = [Vn;
    zeros(size(s.F,1),s.nX),-ones(size(s.F,1),s.nAlpha),zeros(size(s.F,1),s.nBeta),s.F];
v = [vn;ones(size(s.F,1),1)];

[V,v] = rowReduce(V,v);

fprintf('Projecting %d inequalities: ', length(v));
[Aout, bout] = polyProject(V,v,s.nU);
fprintf('Projection has %d inequalities.\n',length(bout));