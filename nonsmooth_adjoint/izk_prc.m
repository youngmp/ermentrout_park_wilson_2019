tstart = 0;
tfinal = 500;
tmaster = [];
ymaster = [];
y0 = [-65;.492421];

a = 0.02;
b = 0.2;
d = 8;
refine=1;

options = odeset('Events',@SpikeEvent,'RelTol',1e-8,'AbsTol',1e-10);

tout = tstart;
yout = y0.';
teout = [];
yeout = [];
ieout = [];

time = 0;
firstflag = 0;
while time < tfinal
    % Solve until the first terminal event.
    [t,y,te,ye,ie] = ode23(@(t,y) izkF(t,y,a,b),[tstart tfinal],y0,options);
  
    % Accumulate output.  This could be passed out as output arguments.
    nt = length(t);
    tout = [tout; t(2:nt)];
    
    yout = [yout; y(2:nt,:)];
    teout = [teout; te];          % Events at tstart are never reported. 
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    
    % save the limit cycle solution
    if firstflag == 0
        lc_y = y;
        lc_t = t;
        
        firstflag = 1;
    end
    
    % Set the new initial conditions with v,u reset.
    y0(1) = -65;
    y0(2) = y(end,2) + d;
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation.  'refine' is 4 by default.
    options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
                             'MaxStep',t(nt)-t(1));
    tstart = t(nt);
    
    time = tstart;
    tmaster = [tmaster;t];
    ymaster = [ymaster;y];
    
end

% phases to perturb
phases = linspace(0,lc_t(end)-.01,100);
prc = zeros(length(phases),1);

% set initial condition with this idx
dx = 0;%1e-2;
dy = 1e-2;

if (dy ~= 0) && (dx ~= 0)
    norm = dy+dx;
    fname_suffix = 'dydx';
elseif dx ~= 0
    norm = dx;
    fname_suffix = 'dx';
elseif dy ~= 0
    norm = dy;
    fname_suffix = 'dy';
end

for i=1:length(prc)

    % choose phase at which to apply x or y pert
    phase = phases(i);
    [C,min_idx] = min(abs(lc_t-phase)); % find index of phase

    

    y02 = lc_y(min_idx,:);

    y02(1,:) = y02(1,:) + [dx dy];

    tstart2 = lc_t(min_idx);

    tmaster2 = [];
    ymaster2 = [];

    options = odeset('Events',@SpikeEvent,'RelTol',1e-8,'AbsTol',1e-10);

    tout2 = tstart2;
    yout2 = y02;
    teout2 = [];
    yeout2 = [];
    ieout2 = [];


    time2 = tstart2;
    while time2 < tfinal
        % Solve until the first terminal event.
        [t2,y2,te2,ye2,ie2] = ode23(@(t,y) izkF(t,y,a,b),[tstart2 tfinal],y02,options);

        % Accumulate output.  This could be passed out as output arguments.
        nt2 = length(t2);
        tout2 = [tout2; t2(2:nt2)];

        yout2 = [yout2; y2(2:nt2,:)];
        teout2 = [teout2; te2];          % Events at tstart are never reported. 
        yeout2 = [yeout2; ye2];
        ieout2 = [ieout2; ie2];

        % Set the new initial conditions with v,u reset.
        y02(1) = -65;
        y02(2) = y2(end,2) + d;
        % A good guess of a valid first timestep is the length of the last valid
        % timestep, so use it for faster computation.  'refine' is 4 by default.
        options = odeset(options,'InitialStep',t2(nt2)-t2(nt2-refine),...
                                 'MaxStep',t2(nt2)-t2(1));
        tstart2 = t2(nt2);

        time2 = tstart2;
        tmaster2 = [tmaster2;t2];
        ymaster2 = [ymaster2;y2];

    end

    %phase difference
    prc(i) = (teout(end)-teout2(end))/norm;

end

figure()
plot(prc)



csvwrite(['matlabprc_' fname_suffix '.csv'],prc)

%{
figure()
plot(tmaster,ymaster(:,1));hold on;
plot(tmaster2,ymaster2(:,1));hold on;

plot(teout,yeout(:,1),'ro'); hold on;
plot(teout2,yeout2(:,1),'ro'); hold off;
%}

%plot(lc_t,lc_y(:,1))

%plot(t,y)