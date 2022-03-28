clear all;
close all;
clc;

% Physical Parameters
global g rho;
g = 9.81;
rho = 1000;

% Create Mesh Object
xL = -1;
xU = 10;
yL = -1;
yU = 1;
nX = 44;
nY = 12;
nVar = 3;
mesh = Mesh(xL, xU, yL, yU, nX, nY, nVar);
% mesh.PlotMesh();

% Define Bathymetry
x = linspace(xL, xU, nX);
y = linspace(yL, yU, nY);
[xMesh, yMesh] = meshgrid(x, y);

% Define Time Domain Information
tstart = 0;
tend = 10;
CFL = 0.4;
nT = 100000;

% Initial Conditions
for i = 1:mesh.nCells
    mesh.data(:,i) = SWEIC(mesh.xCenters(i), mesh.yCenters(i));
end
% mesh.Plot();

% Time Stepping
time = tstart;
fig = figure;
for n = 1:nT
    % Compute time step
    %    Compute eigenvalues of system
    amax = 0;
    for i = 1:mesh.nCells
        amax = max(amax, max(abs(SWELambda(mesh.data(:,i), [0,0]))));
    end
    
    %    Compute dt
    dt = CFL * min(mesh.incircles) / amax;
    
    %    Ensure time stepping ends
    if (time + dt > tend)
        dt = tend - time;
    end
    
    %    Stop criteria
    if (time >= tend)
        break;
    end
    
    % Update solution
    qNew = mesh.data;
    for i = 1:mesh.nCells
        for iEdge = 1:3
            j = mesh.neighbors(i, iEdge);
            nij = mesh.normals(:, iEdge, i);
            if (j < 0)
                qi = mesh.data(:,i);
                qj = qi;
                qj(2:3) = qj(2:3) - 2*(nij(1)*qi(2) + nij(2)*qi(3))*nij;
            else
                qi = mesh.data(:,i);
                qj = mesh.data(:,j);
            end
            smax = max(max(abs(SWELambda(qi, nij))), max(abs(SWELambda(qj, nij))));
            Fn = 0.5*(SWEFlux(qi) + SWEFlux(qj))*nij - 0.5*smax*(qj - qi);
            qNew(:,i) = qNew(:,i) - dt/mesh.areas(i)*mesh.edgeLengths(iEdge,i)*Fn;
        end
        
        % Source Term
        qNew(:,i) = qNew(:,i) + (dt/2)*SWESource(mesh.xCenters(i), mesh.yCenters(i), time, qi);
    end
    
    % Update time and solution
    time = time + dt;
    mesh.data = qNew;
    
    % Integrate over domain for verification
    volume = (xU - xL) * (yU - yL) * sum(mesh.data(1,:));
    mass = volume * rho;
    
    % Plot solution
    if (mod(n, 5) == 0) || (n == 1)
        nodeData = mesh.CenterToNode();
        s = trisurf(mesh.tri, mesh.xGrid, mesh.yGrid, nodeData(1,:));
        set(s, 'EdgeColor', [0 0 0]);
        set(s, 'EdgeAlpha', 0.2);
        set(s, 'FaceColor', 'interp');
        colormap turbo;
        title(sprintf('Time = %1.2f [s], Mass = %4.2f [L]', time, volume / 1000));
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Height [m]');
        view(9, 16.2);
        axis([xL xU yL yU 0 3]);
        drawnow

%         if (mod(n, 50) == 0) || (n == 1)
%             fprintf("Outputing frame %i\n", n);
%             filename = sprintf("frame%i.png", n);
%             saveas(gcf, filename);
%         end
    end
    
    % Create animation
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if n == 1
        imwrite(imind, cm, 'animation.gif', 'gif', 'Loopcount', inf);
    else
        imwrite(imind, cm, 'animation.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
    
end
