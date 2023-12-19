function [me, meLump] = massMatrixBeam(X,rho,A)
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
rotationMatrix = [cosx sinx 0 0 0 0; -sinx cosx 0 0 0 0; 0 0 1 0 0 0; 0 0 0 cosx sinx 0; 0 0 0 -sinx cosx 0; 0 0 0 0 0 1];

%local consistent mass matrix
mLocal = zeros(6,6);
mLocal(1,:) = [rho*A*L/3 0 0 rho*A*L/6 0 0];
mLocal(2,:) = [0 (13/35)*rho*A*L (11/210)*rho*A*L^2 0 (9/70)*rho*A*L -(13/420)*rho*A*L^2];
mLocal(3,:) = [0 (11/210)*rho*A*L^2 (1/105)*rho*A*L^3 0 (13/420)*rho*A*L^2 -(1/140)*rho*A*L^3];
mLocal(4,:) = [rho*A*L/6 0 0 rho*A*L/3 0 0];
mLocal(5,:) = [0 (9/70)*rho*A*L (13/420)*rho*A*L^2 0 (13/35)*rho*A*L -(11/210)*rho*A*L^2];
mLocal(6,:) = [0 -(13/420)*rho*A*L^2 -(1/140)*rho*A*L^3 0 -(11/210)*rho*A*L^2 (1/105)*rho*A*L^3];

%local lump mass matrix
mLocalLump = [(rho*A*L/2) 0 0 0 0 0; 0 (rho*A*L/2) 0 0 0 0; 0 0 (rho*A*L^3/420) 0 0 0;0 0 0 (rho*A*L/2) 0 0; 0 0 0 0 (rho*A*L/2) 0;0 0 0 0 0 (rho*A*L^3/420)];

%global stiffness matrix for the element (consistent and lump)
me = transpose(rotationMatrix)*mLocal*rotationMatrix;
meLump = transpose(rotationMatrix)*mLocalLump*rotationMatrix;

end
