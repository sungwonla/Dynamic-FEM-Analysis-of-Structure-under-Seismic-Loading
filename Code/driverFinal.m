% FINAL PROJECT - Template code for Part 1 - STATIC part

clc;
clf;
close all;
clear all;

% Specify Number of elements
% 12 : Suspension cables
% 13 : Deck
% 4 : Pylons
nel_truss = 12;
nel_beam = 17;
nel = nel_beam + nel_truss;

% Pretension values for the cable elements
f0 = zeros(nel_truss,1);
f0(1) = 7528.46;
f0(2) = 7522.12;
f0(3) = 5896.61;
f0(4) = 5896.61;
f0(5) = 7522.12;
f0(6) = 7528.46;
f0(7) = 7528.46;
f0(8) = 7522.12;
f0(9) = 5896.61;
f0(10) = 7522.12;
f0(11) = 7528.46;
f0(12) = 7528.46;
f0 = -f0;

% Transverse dist. load for the beam elements
w = zeros(nel_beam,1);
w(1:13) = -127.4*1000.0;

%x1,y1,x2, and y2 coordinates of each element are stored in array X
%X(e,1) = x1 of element e,
%X(e,2) = y1 of element e
%X(e,3) = x2.....
X_truss = zeros(nel_truss,4);
X_beam = zeros(nel_beam,4);

% Specify spatial coordinates of truss elements
X_truss(1,:) = [0 30 93 81]; %A-O
X_truss(2,:) = [31 30 93 81]; %B-O
X_truss(3,:) = [62 30 93 81]; %C-O
X_truss(4,:) = [124 30 93 81]; %E-O
X_truss(5,:) = [155 30 93 81]; %F-O
X_truss(6,:) = [186 30 93 81]; %G-O
X_truss(7,:) = [217 30 310 81]; %H-P
X_truss(8,:) = [248 30 310 81]; %I-P
X_truss(9,:) = [279 30 310 81]; %J-P
X_truss(10,:) = [341 30 310 81]; %L-P
X_truss(11,:) = [372 30 310 81]; %M-P
X_truss(12,:) = [403 30 310 81]; %N-P

% Specify spatial coordinates of beam elements
X_beam(1,:) = [0 30 31 30]; %A-B
X_beam(2,:) = [31 30 62 30]; %B-C
X_beam(3,:) = [62 30 93 30]; %C-D
X_beam(4,:) = [93 30 124 30]; %D-E
X_beam(5,:) = [124 30 155 30]; %E-F
X_beam(6,:) = [155 30 186 30]; %F-G
X_beam(7,:) = [186 30 217 30]; %G-H
X_beam(8,:) = [217 30 248 30]; %H-I
X_beam(9,:) = [248 30 279 30]; %I-J
X_beam(10,:) = [279 30 310 30]; %J-K
X_beam(11,:) = [310 30 341 30]; %K-L
X_beam(12,:) = [341 30 372 30]; %L-M
X_beam(13,:) = [372 30 403 30]; %M-N
X_beam(14,:) = [93 0 93 30]; %Q-D
X_beam(15,:) = [93 30 93 81]; %D-O
X_beam(16,:) = [310 0 310 30]; %R-K
X_beam(17,:) = [310 30 310 81]; %K-P

% Young's Modulus : Note that Young's modulus values vary
E_truss = ones(nel_truss,1)*200000*10^6;
E_beam = ones(nel_beam,1)*36500*10^6; %MPa

% Cross sectional area : Note that cross-sections vary
A_truss = zeros(nel_truss,1);
A_beam= zeros(nel_beam,1);

A_truss(1,1) = .1094; %A-O m^2
A_truss(2,1) = .0547; %B-O
A_truss(3,1) = .0547; %C-O
A_truss(4,1) = .0547; %E-O
A_truss(5,1) = .0547; %F-O
A_truss(6,1) = .1094; %G-O
A_truss(7,1) = .1094; %H-P
A_truss(8,1) = .0547; %I-P
A_truss(9,1) = .0547; %J-P
A_truss(10,1) = .0547; %L-P
A_truss(11,1) = .0547; %M-P
A_truss(12,1) = .1094; %N-P

A_beam(1:13,1) = 5.20; %bridge deck
A_beam(14:17,1) = 27.3397; %pylon

% 2nd moment of area
I_beam = zeros(nel_beam,1);
I_beam(1:13,1) = 0.49; %m^4 %bridge deck
I_beam(14:17,1) = 41.4367; %pylon

% Visualize your structure before for unknowns:
% if flagVisualizeInitial is 1 then VISUALIZE
% if flagVisualizeInitial is 0 then skip visualization
flagVisualizeUndeformed = 1;
if flagVisualizeUndeformed == 1
    figure(1)
    for ele=1:nel_truss
        % Undeformed
        p1 = [X_truss(ele,1) X_truss(ele,2)];           % First Point
        p2 = [X_truss(ele,3) X_truss(ele,4)];           % Second Point
        dp = p2-p1;                         % Difference
        q = quiver(p1(1),p1(2),dp(1),dp(2),0);
        q.LineStyle = '-';
        q.Color = 'blue';
        q.AlignVertexCenters = 'on';
        grid
        text(p1(1),p1(2), sprintf('Node:(%.0f,%.0f)',p1))
        text(p2(1),p2(2), sprintf('Node:(%.0f,%.0f)',p2))
        text((p1(1)+p2(1))*0.5,(p1(2)+p2(2))*0.5, sprintf('Truss%.0f ',ele))
        hold on
    end
    for ele=1:nel_beam
        % Undeformed
        p1 = [X_beam(ele,1) X_beam(ele,2)];           % First Point
        p2 = [X_beam(ele,3) X_beam(ele,4)];           % Second Point
        dp = p2-p1;                         % Difference
        q = quiver(p1(1),p1(2),dp(1),dp(2),0);
        q.LineStyle = '-';
        q.Color = 'blue';
        q.AlignVertexCenters = 'on';
        grid
        text(p1(1),p1(2), sprintf('Node:(%.0f,%.0f)',p1))
        text(p2(1),p2(2), sprintf('Node:(%.0f,%.0f)',p2))
        text((p1(1)+p2(1))*0.5,(p1(2)+p2(2))*0.5, sprintf('Beam:%.0f ',ele))
        hold on
    end
    legend('Undeformed Structure');
    xlabel('Global X');
    ylabel('Global Y');
end


%Total number of DOF
ndof = 18*3;

%Total number of free DOF
nfdof = ndof-8;

%Total number of restrained DOF
nddof = ndof - nfdof;

%Specify number of DOFs for each type
edof_beam = 6;
edof_truss = 4;


% ID array. Maps local element DOF to global DOF. Make sure that you put the
% free DOF first and the restrained ones at the end in your own DOF labelling.
% Also, note that ID_beam is of size [nel_beam x edof_beam =6].
% In order to visualize deform structure, please place rotational DOFs to 3rd and 6th entries.
ID_truss = ones(nel_truss,edof_truss);
ID_beam = ones(nel_beam,edof_beam);

%Truss members
ID_truss(1,:) = [1 53 41 42]; %A-O
ID_truss(2,:) = [3 4 41 42]; %B-O
ID_truss(3,:) = [6 7 41 42]; %C-O
ID_truss(4,:) = [12 13 41 42]; %E-O
ID_truss(5,:) = [15 16 41 42]; %F-O
ID_truss(6,:) = [18 19 41 42]; %G-O
ID_truss(7,:) = [21 22 44 45]; %H-P
ID_truss(8,:) = [24 25 44 45]; %I-P
ID_truss(9,:) = [27 28 44 45]; %J-P
ID_truss(10,:) = [33 34 44 45]; %L-P
ID_truss(11,:) = [36 37 44 45]; %M-P
ID_truss(12,:) = [39 54 44 45]; %N-P

% Beams including deck elements and pylons
ID_beam(1,:) = [1 53 2 3 4 5]; %A-B
ID_beam(2,:) = [3 4 5 6 7 8]; %B-C
ID_beam(3,:) = [6 7 8 9 10 11]; %C-D
ID_beam(4,:) = [9 10 11 12 13 14]; %D-E
ID_beam(5,:) = [12 13 14 15 16 17]; %E-F
ID_beam(6,:) = [15 16 17 18 19 20]; %F-G
ID_beam(7,:) = [18 19 20 21 22 23]; %G-H
ID_beam(8,:) = [21 22 23 24 25 26]; %H-I
ID_beam(9,:) = [24 25 26 27 28 29]; %I-J
ID_beam(10,:) = [27 28 29 30 31 32]; %J-K
ID_beam(11,:) = [30 31 32 33 34 35]; %K-L
ID_beam(12,:) = [33 34 35 36 37 38]; %L-M
ID_beam(13,:) = [36 37 38 39 54 40]; %M-N
ID_beam(14,:) = [47 48 49 9 10 11]; %Q-D
ID_beam(15,:) = [9 10 11 41 42 43]; %D-O
ID_beam(16,:) = [50 51 52 30 31 32]; %R-K
ID_beam(17,:) = [30 31 32 44 45 46]; %K-P


%Initialize stiffness matrix K
% Global Stiffness Matrix
K = zeros(ndof,ndof);

%Initialize force vector P (external nodal forces) and Pf0 (fixed end forces on nodes)
%Put you code here to populate array P based on the problem statement
%You can only populate up to nfdof index. The indices nfdof+1:ndof
%correspond to the restrained DOF and are in fact the reaction forces that
%will be calculated later on
P = zeros(ndof,1);
Pf0 = zeros(ndof,1);


%Populating global stiffness matrix
% Assemble truss elements
for e=1:nel_truss
    [ke_truss, fe_truss] = localStiffnessTruss(X_truss(e,:),E_truss(e,1),A_truss(e,1),f0(e,1));
    for i=1:edof_truss  %Loopoing over local dof
        for j=1:edof_truss  %Looping over local dof
            K(ID_truss(e,i),ID_truss(e,j)) = K(ID_truss(e,i),ID_truss(e,j)) + ke_truss(i,j);
        end
        Pf0(ID_truss(e,i),1) = Pf0(ID_truss(e,i),1) + fe_truss(i,1);
    end
end

% Assemble beam elements
for e=1:nel_beam
    [ke_beam, fe_beam]= localStiffnessBeam(X_beam(e,:),E_beam(e,1),I_beam(e,1),A_beam(e,1),w(e,1));
    for i=1:edof_beam  %Loopoing over local dof
        for j=1:edof_beam  %Looping over local dof
            K(ID_beam(e,i),ID_beam(e,j)) = K(ID_beam(e,i),ID_beam(e,j)) + ke_beam(i,j);
        end
        Pf0(ID_beam(e,i),1) = Pf0(ID_beam(e,i),1) + fe_beam(i,1);
    end
end


%We can split the stiffness matrix into four submatrices. Kff, Kfd, Kdf, Kdd
%Kff has dimensions nfdof x nfdof
%Kfd has dimensions nddof x nfdof
%Kdf has dimensions nfdof x nddof
%Kdd has dimensions nddof x nddof
%Extracting the submatrices from matrix K
Kff = K(1:nfdof,1:nfdof);
Kfd = K(1:nfdof,nfdof+1:ndof);
Kdf = K(nfdof+1:ndof,1:nfdof);
Kdd = K(nfdof+1:ndof,nfdof+1:ndof);

%Extracting Pf from P
Pf = P(1:nfdof,1);

%Extracting Pf0f and Pf0d from Pf0 (fixed end forces on structural level)
Pf0f = Pf0(1:nfdof,1);
Pf0d = Pf0(nfdof+1:ndof,1);

%Initialize array of displacements u
u = zeros(ndof,1);

%Extracting ud from u
ud = u(nfdof+1:ndof,1);

%Solving for the displacements uf
uf = Kff \ (Pf - Pf0f - Kfd*ud);

%Solving for the reaction forces
Pd = Kdf*uf + Kdd*ud + Pf0d;


% Visualize your deformed structure:
% if flagVisualizeInitial is 1 then VISUALIZE
% if flagVisualizeInitial is 0 then skip visualization
flagVisualizeDeformed = 1;
% Manually adjust below amplification factor so that you could visualize your results
amplificationFactor = 50.0;
if flagVisualizeDeformed == 1
    % Retrieve full u
    u = [uf(:)' ud(:)'];
    % Define a 2D arrays containing updated positions of the nodes
    X_def_truss = zeros(nel,4);
    X_def_vis_truss = zeros(nel,4);

    X_def_beam = zeros(nel,4);
    X_def_vis_beam = zeros(nel,4);
    % Loop over each dof to update positions
    for ele =1:nel_truss
        for i=1:edof_truss
            X_def_truss(ele,i) = X_truss(ele,i) + u(ID_truss(ele,i));
            X_def_vis_truss(ele,i) = X_truss(ele,i) + amplificationFactor*u(ID_truss(ele,i));
        end
    end
    % Only take into account displacements.
    for ele =1:nel_beam
        X_def_beam(ele,1) = X_beam(ele,1) + u(ID_beam(ele,1));
        X_def_vis_beam(ele,1) = X_beam(ele,1) + amplificationFactor*u(ID_beam(ele,1));

        X_def_beam(ele,2) = X_beam(ele,2) + u(ID_beam(ele,2));
        X_def_vis_beam(ele,2) = X_beam(ele,2) + amplificationFactor*u(ID_beam(ele,2));

        X_def_beam(ele,3) = X_beam(ele,3) + u(ID_beam(ele,5));
        X_def_vis_beam(ele,3) = X_beam(ele,3) + amplificationFactor*u(ID_beam(ele,4));

        X_def_beam(ele,4) = X_beam(ele,4) + u(ID_beam(ele,6));
        X_def_vis_beam(ele,4) = X_beam(ele,4) + amplificationFactor*u(ID_beam(ele,5));
    end
    figure(2)
    for ele=1:nel_truss
        hold on
        % Undeformed
        p1 = plot([X_truss(ele,1) X_truss(ele,3)], [X_truss(ele,2) X_truss(ele,4)],'--b');
        % Deformed
        p2 = plot([X_def_vis_truss(ele,1) X_def_vis_truss(ele,3)], [X_def_vis_truss(ele,2) X_def_vis_truss(ele,4)]...
            , 'k','LineWidth',1.5);
        legend([p1 p2],'Undeformed','Deformed');
    end
    for ele=1:nel_beam
        hold on
        % Undeformed
        p1 = plot([X_beam(ele,1) X_beam(ele,3)], [X_beam(ele,2) X_beam(ele,4)],'--r');
        % Deformed
        p2 = plot([X_def_vis_beam(ele,1) X_def_vis_beam(ele,3)], [X_def_vis_beam(ele,2) X_def_vis_beam(ele,4)]...
            , 'k','LineWidth',1.5);
        legend([p1 p2],'Undeformed','Deformed');
        hold off
    end
    title(['Deformed Structure: Amplification factor: ',num2str(amplificationFactor)]);
    xlabel('Global X');
    ylabel('Global Y');
end



%Part 2: Dynamics
rho_truss = ones(nel_truss,1)*8050; %kg/m^3
rho_beam = ones(nel_beam,1)*2400; %kg/m^3

%Populating global mass matrix (consistent and lump)
M = zeros(ndof,ndof);
MLump = zeros(ndof,ndof);
% Assemble truss elements
for e=1:nel_truss
    [mass_truss,mass_truss_lump] = massMatrixTruss(X_truss(e,:),rho_truss(e,1),A_truss(e,1));
    for i=1:edof_truss  %Loopoing over local dof
        for j=1:edof_truss  %Looping over local dof
            M(ID_truss(e,i),ID_truss(e,j)) = M(ID_truss(e,i),ID_truss(e,j)) + mass_truss(i,j);
            MLump(ID_truss(e,i),ID_truss(e,j)) = MLump(ID_truss(e,i),ID_truss(e,j)) + mass_truss_lump(i,j);
        end
    end
end

% Assemble beam elements
for e=1:nel_beam
    [mass_beam,mass_beam_lump] = massMatrixBeam(X_beam(e,:),rho_beam(e,1),A_beam(e,1));
    for i=1:edof_beam  %Loopoing over local dof
        for j=1:edof_beam  %Looping over local dof
            M(ID_beam(e,i),ID_beam(e,j)) = M(ID_beam(e,i),ID_beam(e,j)) + mass_beam(i,j);
            MLump(ID_beam(e,i),ID_beam(e,j)) = MLump(ID_beam(e,i),ID_beam(e,j)) + mass_beam_lump(i,j);
        end
    end
end

%We can split the mass matrix into four submatrices. Mff, Mfd, Mdf, Mdd
%Mff has dimensions nfdof x nfdof
%Mfd has dimensions nddof x nfdof
%Mdf has dimensions nfdof x nddof
%Mdd has dimensions nddof x nddof
%Extracting the submatrices from matrix M
Mff = M(1:nfdof,1:nfdof);
Mfd = M(1:nfdof,nfdof+1:ndof);
Mdf = M(nfdof+1:ndof,1:nfdof);
Mdd = M(nfdof+1:ndof,nfdof+1:ndof);
%Extracting the submatrices from matrix MLump
MffLump = MLump(1:nfdof,1:nfdof);
MfdLump = MLump(1:nfdof,nfdof+1:ndof);
MdfLump = MLump(nfdof+1:ndof,1:nfdof);
MddLump = MLump(nfdof+1:ndof,nfdof+1:ndof);

%forcing vector from distributed loads
forcingVector = Pf - Pf0f;

%defining numerical approximation parameters
rho_inf = 0.5;
am = (2-rho_inf)/(1+rho_inf);
af = 1/(1+rho_inf);
gamma = 0.5+am-af;
beta = 0.25*(1+am-af)^2;

%time step
dt = 0.02;
tfinal = 40;
t = 0:dt:tfinal;

%initial conditions, from static part 1
%consistent mass matrix
ufAcel = zeros(nfdof,1);
ufVel = zeros(nfdof,1);
ufDisp = uf;
%lump mass matrix
ufAcelLump = zeros(nfdof,1);
ufVelLump = zeros(nfdof,1);
ufDispLump = uf;

%prescribed DOFs, representing earthquake movement
udAcel = readmatrix("acc.txt");
udVel = readmatrix("vel.txt");
udDisp = readmatrix("disp.txt");

%u vector to plot (consistent and lump mass matrix)
uConsistent = u;
uLump = u;

%matrix storing reaction forces at prescribed DOFs
PdNew = zeros(nddof,length(t));
PdNewLump = zeros(nddof,length(t));

%matrix storing acceleration components at nodes D and K
accelD = zeros(1,length(t));
accelDLump = zeros(1,length(t));
accelK = zeros(1,length(t));
accelKLump = zeros(1,length(t));

for i = 1:length(t)-1
    %consistent mass matrix
    %predictor stage
    ufAcelNextPredict = ((gamma-1)/gamma)*ufAcel;
    ufVelNextPredict = ufVel;
    ufDispNextPredict = ufDisp + dt*ufVel + (dt^2/2)*((1-2*beta)*ufAcel + (2*beta)*ufAcelNextPredict);
    %corrector stage
    ufAcelAlpha = ufAcel + am*(ufAcelNextPredict-ufAcel);
    ufVelAlpha = ufVel + af*(ufVelNextPredict-ufVel);
    ufDispAlpha = ufDisp + af*(ufDispNextPredict-ufDisp);
    %ualpha at prescribed DOFs, from given data
    udAcelAlpha = udAcel(i) + am*(udAcel(i+1)-udAcel(i));
    udVelAlpha = udVel(i) + af*(udVel(i+1)-udVel(i));
    udDispAlpha = udDisp(i) + af*(udDisp(i+1)-udDisp(i));
    udAcelAlphaMatrix = [udAcelAlpha 0 0 udAcelAlpha 0 0 0 0]';
    udVelAlphaMatrix = [udVelAlpha 0 0 udVelAlpha 0 0 0 0]';
    udDispAlphaMatrix = [udDispAlpha 0 0 udDispAlpha 0 0 0 0]';

    %lump mass matrix
    %predictor stage
    ufAcelNextPredictLump = ((gamma-1)/gamma)*ufAcelLump;
    ufVelNextPredictLump = ufVelLump;
    ufDispNextPredictLump = ufDispLump + dt*ufVelLump + (dt^2/2)*((1-2*beta)*ufAcelLump + (2*beta)*ufAcelNextPredictLump);
    %corrector stage
    ufAcelAlphaLump = ufAcelLump + am*(ufAcelNextPredictLump-ufAcelLump);
    ufVelAlphaLump = ufVelLump + af*(ufVelNextPredictLump-ufVelLump);
    ufDispAlphaLump = ufDispLump + af*(ufDispNextPredictLump-ufDispLump);

    %Animation for Consistent Mass Matrix
    % Call the following routine after corrector step in every time step
    % Create an animation
    flagVisualizeAnimation = 1;
    figure_object = figure(3); % Create figure for animation
    filename = 'animationEQConsistent.gif';
    amplificationFactor = 5;
    if flagVisualizeAnimation == 1
        % Append to animation
        visualize(uConsistent,nel_truss, nel_beam, ID_truss, ID_beam, X_truss, X_beam, figure_object, amplificationFactor);
        drawnow
        pause(0.000)
        % Capture the plot as an image
        frame = getframe(figure_object);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if i == 1 %t == 0.0
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append');
        end
    end
    %Animation for Lump Mass Matrix
    % Call the following routine after corrector step in every time step
    % Create an animation
    flagVisualizeAnimation = 1;
    figure_object = figure(4); % Create figure for animation
    filename = 'animationEQLump.gif';
    amplificationFactor = 5;
    if flagVisualizeAnimation == 1
        % Append to animation
        visualize(uLump,nel_truss, nel_beam, ID_truss, ID_beam, X_truss, X_beam, figure_object, amplificationFactor);
        drawnow
        pause(0.000)
        % Capture the plot as an image
        frame = getframe(figure_object);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if i == 1 %t == 0.0
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append');
        end
    end
    % TIME INTEGRATION LOOP ENDS

    %collocation, solve for acceleration increment (consistent and lump mass matrix)
    accelIncrement = (am*Mff + af*beta*(dt^2)*Kff)\(forcingVector - (Mff*ufAcelAlpha+Kff*ufDispAlpha) - (Mfd*udAcelAlphaMatrix+Kfd*udDispAlphaMatrix));
    accelIncrementLump = (am*MffLump + af*beta*(dt^2)*Kff)\(forcingVector - (MffLump*ufAcelAlphaLump+Kff*ufDispAlphaLump) - (MfdLump*udAcelAlphaMatrix+Kfd*udDispAlphaMatrix));

    %update solution
    %consistent mass matrix
    ufAcel = ufAcelNextPredict + accelIncrement;
    ufVel = ufVelNextPredict + gamma*dt*accelIncrement;
    ufDisp = ufDispNextPredict + beta*(dt^2)*accelIncrement;
    %lump mass matrix
    ufAcelLump = ufAcelNextPredictLump + accelIncrementLump;
    ufVelLump = ufVelNextPredictLump + gamma*dt*accelIncrementLump;
    ufDispLump = ufDispNextPredictLump + beta*(dt^2)*accelIncrementLump;

    %post-processing to get reaction forces and accelerations at D and K,
    %consistent and lump mass matrix
    udAcelMatrix = [udAcel(i+1) 0 0 udAcel(i+1) 0 0 0 0]';
    udDispMatrix = [udDisp(i+1) 0 0 udDisp(i+1) 0 0 0 0]';
    PdNew(:,i) = Pf0d + Mdf*ufAcel + Mdd*udAcelMatrix + Kdf*ufDisp + Kdd*udDispMatrix;
    accelD(i+1) = ufAcel(9);
    accelK(i+1) = ufAcel(30);
    PdNewLump(:,i) = Pf0d + MdfLump*ufAcelLump + MddLump*udAcelMatrix + Kdf*ufDispLump + Kdd*udDispMatrix;
    accelDLump(i+1) = ufAcelLump(9);
    accelKLump(i+1) = ufAcelLump(30);

    uConsistent = [ufDisp;udDispMatrix];
    uLump = [ufDispLump;udDispMatrix];
    fprintf('t = %e\n', t(i));
end

%%%%UNCOMMENT TO VIEW PLOTS%%%%%%%%%%%

%plotting acceleration components at nodes D and K, consistent and lump mass matrix
%plotting support reactions (1-3 are node Q x-y-moment, 4-6 are node R x-y-moment, 7 is node A y direction, 8
%is node N y direction)
figure(5);
tiledlayout(2,2)
ax1 = nexttile;
plot(t,accelD);
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
title('X Acceleration at Node D, Consistent Mass Matrix');
ax2 = nexttile;
plot(t,accelK);
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
title('X Acceleration at Node K, Consistent Mass Matrix');
ax3 = nexttile;
plot(t,accelDLump);
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
title('X Acceleration at Node D, Lumped Mass Matrix');
ax4 = nexttile;
plot(t,accelKLump);
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
title('X Acceleration at Node K, Lumped Mass Matrix');
sgtitle('X Acceleration Components at Nodes D and K vs. Time');
linkaxes([ax1 ax2 ax3 ax4],'xy')
ax1.XLim = [0 tfinal];

figure(6);
tiledlayout(4,2)
ax5 = nexttile;
plot(t,PdNew(1,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (X) at node Q');
ax6 = nexttile;
plot(t,PdNew(2,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (Y) at node Q');
ax7 = nexttile;
plot(t,PdNew(3,:));
xlabel('Time (s)');
ylabel('M (N-m)');
title('Moment at node Q');
ax8 = nexttile;
plot(t,PdNew(4,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (X) at node R');
ax9 = nexttile;
plot(t,PdNew(5,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (Y) at node R');
ax10 = nexttile;
plot(t,PdNew(6,:));
xlabel('Time (s)');
ylabel('M (N-m)');
title('Reaction Moment at node R');
ax11 = nexttile;
plot(t,PdNew(7,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (Y) at node A');
ax12 = nexttile;
plot(t,PdNew(8,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (Y) at node N');
sgtitle('Support Reaction Forces vs. Time, Consistent Mass Matrix');
linkaxes([ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12],'xy')
ax5.XLim = [0 tfinal];

figure(7);
tiledlayout(2,4)
ax13 = nexttile;
plot(t,PdNewLump(1,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (X) at node Q');
ax14 = nexttile;
plot(t,PdNewLump(2,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (Y) at node Q');
ax15 = nexttile;
plot(t,PdNewLump(3,:));
xlabel('Time (s)');
ylabel('M (N-m)');
title('Reaction Moment at node Q');
ax16 = nexttile;
plot(t,PdNewLump(4,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (X) at node R');
ax17 = nexttile;
plot(t,PdNewLump(5,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (Y) at node R');
ax18 = nexttile;
plot(t,PdNewLump(6,:));
xlabel('Time (s)');
ylabel('M (N-m)');
title('Reaction Moment at node R');
ax19 = nexttile;
plot(t,PdNewLump(7,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (Y) at node A');
ax20 = nexttile;
plot(t,PdNewLump(8,:));
xlabel('Time (s)');
ylabel('F (N)');
title('Force (Y) at node N');
sgtitle('Support Reaction Forces vs. Time, Lumped Mass Matrix');
linkaxes([ax13 ax14 ax15 ax16 ax17 ax18 ax19 ax20],'xy')
ax13.XLim = [0 tfinal];
