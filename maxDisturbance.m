% This script computes the maximal disturbance the system can cope with.

clear all
close all 
load Hor3Const

alphaMax = 2;
r = 2;

[OBJ,RHS] = rowReduce(Lam{1}(:,[1,2,4]),lam{1} - Lam{1}(:,3)*alphaMax);

model.obj = [0,0,1];
model.A = sparse(Lam{1}(:,[1,2,4]));
model.sense = '<';
model.rhs = lam{1} - Lam{1}(:,3)*alphaMax;
model.modelsense = 'max';
model.varnames = {'x1', 'x2', 'beta'};

param.OutputFlag = 0;

res = gurobi(model,param);

clc;
fprintf('For alpha = %.3f the maximal beta is %.3f. \n',...
        alphaMax, res.objval);

clear model;
    
    
model.quadcon(1).Qc = sparse(eye(2));
model.quadcon(1).q = [0,0];
model.quadcon(1).rhs = r^2;
model.A = sparse(zeros(1,2));
model.sense = '<';
model.rhs = 1;

result = zeros(size(RHS));

for i = 1:length(RHS)
    model.obj = -OBJ(i,1:2);
    res = gurobi(model,param);
    if and(OBJ(i,3)==0,res.objval>RHS(i))
        result(i) = -inf;
    elseif OBJ(i,3)<=0
        result(i) = inf;
    else
        result(i) = (RHS(i)+res.objval)/OBJ(i,3);
    end
end

betaMax = min(result);

if betaMax>0
fprintf('For alpha = %.3f the maximal beta that contains the disk of radius %.3f is %.3f. \n',...
        alphaMax, r, betaMax);
else
    fprintf('The disk of radius %.3f is not contained in the invariant set for alpha = %.3f.\n',...
        r,alphaMax)
end

[V,v] = rowReduce(Lam{1}(:,1:2),lam{1}-Lam{1}(:,3)*alphaMax-Lam{1}(:,4)*betaMax);
figure(1)
P = Polyhedron(V,v);
plot(P)
hold on
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
N = 300;

[X,Y] = meshgrid(linspace(xlim(1),xlim(2),N),linspace(ylim(1),ylim(2),N));
Z = zeros(N);
for i = 1:N
    for j = 1:N
        if (X(i,j)^2+Y(i,j)^2)<=r^2
            Z(i,j) = 1;
        else
            Z(i,j) = 0;
        end
    end
end

[C,h] = contourf(X,Y,Z,'Linewidth',0.1);
allH = allchild(h);
valueToHide = 1;
patchValues = cell2mat(get(allH,'UserData'));
patchesToHide = patchValues == valueToHide;
set(allH(patchesToHide),'FaceColor','k','FaceAlpha',0.8);
plot(P,'alpha',.2)
hold off