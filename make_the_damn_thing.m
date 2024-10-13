%This is an implementation of swarmalators as found in
% https://en.wikipedia.org/wiki/Swarmalators

%Set up
J = 1; K = -0.75; %Active phase wave
% J = 1; K = 0; %Static phase wave
% J = 1; K = -0.1; %Splintered phase wave
num_swarmers = 200;
%Settings for gif generation 
num_frames = 200; %For eventual plot
delay_time = 1/10; %in seconds
filename = 'swarmers.gif';

%Make initial conditions
pos0 = rand(num_swarmers,2);
pos0 = stretch_to_interval(pos0,[0 1],[-2 2]);

phase0 = rand(num_swarmers,1);
phase0 = stretch_to_interval(phase0,[0 1],[0 2*pi]);


u0 = [pos0(:,1); pos0(:,2); phase0];

%Solve!
Tfin = 60;
[t,u] = ode45(@(t,u)swarmalator_ode(t,u,num_swarmers,J,K), [0 Tfin], u0);

%Get output ready for easy plotting
times = linspace(0,Tfin,num_frames);
u_interp = interp1(t,u,times);
pos = zeros(num_frames,num_swarmers,2);
phase = zeros(num_frames,num_swarmers);
%I suspect there is a more vector-savvy way to do this
for i = 1:num_swarmers
    swarmer_pos_idx = i + [0,1]*num_swarmers;
    pos(:,i,:) = u_interp(:,swarmer_pos_idx);
    phase(:,i) = u_interp(:,2*num_swarmers+i);
end
phase = stretch_to_interval(phase,[0 2*pi],[0 1]);
phase = mod(phase,1);

%Plot!
figure; hold on
%Get axis limits for whole timeline
XLim = [min(min(pos(:,:,1))) max(max(pos(:,:,1)))];
YLim = [min(min(pos(:,:,2))) max(max(pos(:,:,2)))];
axis([XLim YLim]);
axis off
set(gcf,'Color',0.8*[1 1 1]);
swarmer_ps = zeros(1,num_swarmers);
for j = 1:num_swarmers
    swarmer_ps(j) = scatter(pos(1,j,1), pos(1,j,2), 60, 'ko');
    set(swarmer_ps(j),'MarkerFaceColor',hsv2rgb([phase(1,j) 1 1]));
    set(swarmer_ps(j),'MarkerFaceAlpha',0.5);
end
drawnow
frame = getframe(gcf);
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delay_time);
for i = 1:num_frames
    for j = 1:num_swarmers
        set(swarmer_ps(j),'XData',pos(i,j,1));
        set(swarmer_ps(j),'YData',pos(i,j,2));
        set(swarmer_ps(j),'MarkerFaceColor',hsv2rgb([phase(i,j) 1 1]));
    end
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delay_time);
end

function stretched = stretch_to_interval(u,I1,I2)
%Values are mapped from interval I1 to I2.
    u = (u - I1(1)) / diff(I1);
    stretched = diff(I2)*u + I2(1);
end

function du = swarmalator_ode(~,u,num_swarmers,J,K)
    du = zeros(size(u));

    % Repackage to something fit for human eyes
    pos = zeros(num_swarmers,2);
    for i = 1:num_swarmers
        swarmer_pos_idx = i + [0,1]*num_swarmers;
        pos(i,:) = u(swarmer_pos_idx);
    end
    phase = u(2*num_swarmers+1:end);

    for i = 1:num_swarmers
        %Compare position and phase with neighbors
        others_idx = 1:num_swarmers~=i;
        pos_diffs = pos(others_idx,:) - pos(i,:);
        pos_diffs_norm = vecnorm(pos_diffs,2,2);
        phase_diffs = phase(others_idx,:) - phase(i,:);

        %These are the actual ODE terms
        dpos = 1/num_swarmers * sum( ...
                pos_diffs./pos_diffs_norm.*(1+J*cos(phase_diffs)) ...
                - pos_diffs./pos_diffs_norm.^2, ...
                1 ...
            );
        dphase = K/num_swarmers * sum( ...
                sin(phase_diffs)./pos_diffs_norm, ...
                1 ...
            );

        %Fit into annoying du shape
        du(i) = dpos(1);
        du(i+num_swarmers) = dpos(2);
        du(i+2*num_swarmers) = dphase;
    end
end