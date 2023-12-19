% Routine for generating animations
% u: Displacement vector
% nel_truss: Number of truss elements
% nel_beam: Number of beam elements
% ID array for truss: ID_truss
% ID array for beam:  ID_beam
% Node coordinates of truss: X_truss,
% Node coordinates of beam: X_beam
% Figure object: fig
% Amplification factor for animation: amp
function [] = visualize(u,nel_truss, nel_beam, ID_truss, ID_beam, X_truss, X_beam, fig, amp)
    % Define a 2D arrays containing updated positions of the nodes
    clf(fig);
    X_def_truss = zeros(nel_truss,4);
    X_def_vis_truss = zeros(nel_truss,4);
    edof_truss = 4;
    X_def_beam = zeros(nel_beam,4);
    X_def_vis_beam = zeros(nel_beam,4);
    % Loop over each dof to update positions
    for ele =1:nel_truss
        for i=1:edof_truss
            X_def_truss(ele,i) = X_truss(ele,i) + u(ID_truss(ele,i));
            X_def_vis_truss(ele,i) = X_truss(ele,i) + amp*u(ID_truss(ele,i));
        end
    end
    % Only take into account displacements.
    for ele =1:nel_beam
        X_def_beam(ele,1) = X_beam(ele,1) + u(ID_beam(ele,1));
        X_def_vis_beam(ele,1) = X_beam(ele,1) + amp*u(ID_beam(ele,1));

        X_def_beam(ele,2) = X_beam(ele,2) + u(ID_beam(ele,2));
        X_def_vis_beam(ele,2) = X_beam(ele,2) + amp*u(ID_beam(ele,2));

        X_def_beam(ele,3) = X_beam(ele,3) + u(ID_beam(ele,5));
        X_def_vis_beam(ele,3) = X_beam(ele,3) + amp*u(ID_beam(ele,4));

        X_def_beam(ele,4) = X_beam(ele,4) + u(ID_beam(ele,6));
        X_def_vis_beam(ele,4) = X_beam(ele,4) + amp*u(ID_beam(ele,5));
    end
    %figure(2)
    for ele=1:nel_truss
        % Undeformed
        fig = plot([X_truss(ele,1) X_truss(ele,3)], [X_truss(ele,2) X_truss(ele,4)],'--b');
        axis auto
        hold on
        % Deformed
        fig = plot([X_def_vis_truss(ele,1) X_def_vis_truss(ele,3)], [X_def_vis_truss(ele,2) X_def_vis_truss(ele,4)]...
            , 'k','LineWidth',1.5);
    end
    for ele=1:nel_beam
        fig = plot([X_beam(ele,1) X_beam(ele,3)], [X_beam(ele,2) X_beam(ele,4)],'--r');
        hold on
        % Deformed
        fig = plot([X_def_vis_beam(ele,1) X_def_vis_beam(ele,3)], [X_def_vis_beam(ele,2) X_def_vis_beam(ele,4)]...
            , 'k','LineWidth',1.5);
        
        xlim([-100 500]) 
        ylim([0 150])
    end
    title(['Deformed Structure: Amplification factor: ',num2str(amp)]);
    xlabel('Global X');
    ylabel('Global Y');      
end