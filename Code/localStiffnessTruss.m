function [ke, fe] = localStiffnessTruss(X,E,A,f0)
% Given: spatial coordinates of the nodes X(.,:), Young's Modulus (E),
% cross- sectional area (A) and pre-tensioning value of the member (f0)
% Return: element global stiffness matrix and element fixed end forces

ke = zeros(4,4);

%input variables
deltaX = X(3)-X(1);
deltaY = X(4)-X(2);

%calculate lenghts and cosines/sines
length = sqrt(deltaX^2 + deltaY^2);
cosx = deltaX/length;
sinx = deltaY/length;

%change basis to local/basic systems
rotationMatrix = [cosx sinx 0 0; -sinx cosx 0 0; 0 0 cosx sinx; 0 0 -sinx cosx];
rbmMatrix = [-1 0 1 0];

%global stiffness matrix
ke = transpose(rotationMatrix)*transpose(rbmMatrix)*[E*A/length]*rbmMatrix*rotationMatrix;

%element fixed end forces
felocal = [-f0*length/2 0 -f0*length/2 0]';
fe = transpose(rotationMatrix)*felocal;

end
