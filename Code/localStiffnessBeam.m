function [ke,fe] = localStiffnessBeam(X,E,I,A,w)
% Given: spatial coordinates of the nodes, Young's Modulus, Inertia (I),
% cross- sectional area (A) and distributed load (w)
% Return: element global stiffness matrix and element fixed end forces

%input variables
deltaX = X(3)-X(1);
deltaY = X(4)-X(2);

%calculate lengths and cosines/sines
length = sqrt(deltaX^2 + deltaY^2);
cosx = deltaX/length;
sinx = deltaY/length;

%change basis to local/basic systems
rotationMatrix = [cosx sinx 0 0 0 0; -sinx cosx 0 0 0 0; 0 0 1 0 0 0; 0 0 0 cosx sinx 0; 0 0 0 -sinx cosx 0; 0 0 0 0 0 1];
rbmMatrix = [0 1/length 1 0 -1/length 0; 0 1/length 0 0 -1/length 1; -1 0 0 1 0 0];

kBasic = [4*E*I/length 2*E*I/length 0; 2*E*I/length 4*E*I/length 0; 0 0 E*A/length];
%global stiffness matrix for the element
ke = transpose(rotationMatrix)*transpose(rbmMatrix)*kBasic*rbmMatrix*rotationMatrix;

%element fixed end forces
wx = sinx*w;
wy = cosx*w;
fexlocal = [-0.5*wx*length 0 0 -0.5*wx*length 0 0]';
feylocal = [0 -0.5*wy*length (-1/12)*wy*length^2 0 -0.5*wy*length (1/12)*wy*length^2]';
fexglobal = transpose(rotationMatrix)*fexlocal;
feyglobal = transpose(rotationMatrix)*feylocal;
fe = fexglobal + feyglobal;

end
