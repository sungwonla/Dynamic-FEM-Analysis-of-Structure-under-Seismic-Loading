function [me, meLump] = massMatrixTruss(X,rho,A)
% Given: spatial coordinates of the nodes, density(rho), cross-sectional area (A)
% Return: element global mass matrix (both consistent and lump mass matrix)

%input variables
deltaX = X(3)-X(1);
deltaY = X(4)-X(2);

%calculate Ls and cosines/sines
L = sqrt(deltaX^2 + deltaY^2);
cosx = deltaX/L;
sinx = deltaY/L;

%change basis to local/basic systems
rotationMatrix = [cosx sinx 0 0; -sinx cosx 0 0; 0 0 cosx sinx; 0 0 -sinx cosx];

%local consistent and lump mass matrix
mLocal = rho*A*L*[1/3 0 1/6 0; 0 1/3 0 1/6; 1/6 0 1/3 0; 0 1/6 0 1/3];
mLocalLump = rho*A*L*[1/2 0 0 0; 0 1/2 0 0; 0 0 1/2 0; 0 0 0 1/2];

%global stiffness matrix for the element (consistent and lump)
me = transpose(rotationMatrix)*mLocal*rotationMatrix;
meLump = transpose(rotationMatrix)*mLocalLump*rotationMatrix;

end
