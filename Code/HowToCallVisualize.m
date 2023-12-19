% Call the following routine after corrector step in every time step
% Create an animation
flagVisualizeAnimation = 1;
figure_object = figure(1); % Create figure for animation
filename = 'animationEQ.gif';
amplificationFactor = 5;
if flagVisualizeAnimation == 1
    % Append to animation
    visualize(u,nel_truss, nel_beam, ID_truss, ID_beam, X_truss, X_beam, figure_object, amplificationFactor);
    drawnow
    pause(0.000)
    % Capture the plot as an image 
    frame = getframe(figure_object); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if t == 0.0 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    end 
end
% TIME INTEGRATION LOOP ENDS