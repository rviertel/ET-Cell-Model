function data = ET(ORNtrace, ORNsamplingrate)

    %% set up params
	PARS = parameters;

    ics = [-51.7045    0.0531    0.0604    0.1720    0.0460    0.2084    0.0005    0.1292 0]';

    odeopts = odeset('Events',@spikedetect);

    data = integrator(@vfield_ET, odeopts, ics, PARS, ORNtrace, ORNsamplingrate);
end

function xdot = vfield_ET(t,x,p)
% Vector field for the minimal ET model.
%
%   INPUTS:
%   t -- current time
%   x -- (6,1) vector of current vector values
%   p -- struct containing parameter values in p.ET.param_name format
%
%   OUTPUT:
%   xdot -- derivative w.r.t. time

%% Phase variables

V = x(1);
nK = x(2);
hNaP = x(3);
hH = x(4);
mLVA = x(5);
hLVA = x(6);
wBK = x(7);
Ca = x(8);
nNew = x(9);


% auxiliary quantities for the BK current
theta_wBK = -32 + 59.2*exp(-90*Ca) + 96.7*exp(-470*Ca);
p_wBK = 2.9 + 6.3*exp(-360*Ca);
s = -25.3 + 107.5*exp(-120*Ca);
f = 1/(10*(exp(-(V + 100 - s)/63.6)+exp((-150+(V + 100 - s))/63.6))) - 5.2;



% infinity curves
mNa_inf = 1./(1+exp((V - p.theta_mNa)./p.sigma_mNa));
nK_inf = 1./(1+exp((V - p.theta_nK)./p.sigma_nK));
mNaP_inf = 1./(1+exp((V - p.theta_mNaP)./p.sigma_mNaP));
hNaP_inf = 1./(1+exp((V - p.theta_hNaP)./p.sigma_hNaP));
hH_inf = 1./(1+exp((V - p.theta_hH)./p.sigma_hH));
mLVA_inf = 1/(1+exp((V - p.theta_mLVA)/p.sigma_mLVA));
hLVA_inf = 1/(1+exp((V - p.theta_hLVA)/p.sigma_hLVA));
mHVA_inf = 1/(1+exp((V - p.theta_mHVA)/p.sigma_mHVA));
wBK_inf = 1/(1+exp((V - theta_wBK)/p.sigma_wBK));
mNew_inf = 1./(1+exp((V - p.theta_mNew)/p.sigma_mNew));
nNew_inf = 1./(1+exp((V - p.theta_nNew)/p.sigma_nNew));

% time constants
nNew_tau = 1000./(1.0+exp(-(V+35))) + 1000;
nNew_tau = 1000./(1.0+exp(-(V+35))) + 1000;
nK_tau = p.tau_nK./cosh((V-p.theta_nK)/(2.0*p.sigma_nK));
hNaP_tau = p.tau_hNaP./cosh((V-p.theta_hNaP)/(2.0*p.sigma_hNaP));
hH_tau = p.tau_hH_T*exp(p.delta_hH_T*(V-p.theta_hH_T)/p.sigma_hH_T)/(1+exp((V-p.theta_hH_T)/p.sigma_hH_T));
mLVA_tau = p.tau_mLVA./cosh((V-p.theta_mLVA)/(2.0*p.sigma_mLVA));
hLVA_tau = p.tau_hLVA./cosh((V-p.theta_hLVA)/(2.0*p.sigma_hLVA));
wBK_tau = -(p_wBK - 1.0)*(f - 0.2)/0.8 + p.wBK_base;

% compute values for the currents
INa = p.gNa.*(1-nK)*mNa_inf.^3*(V - p.vNa);
INew = p.gNew*mNew_inf*nNew*(V - p.vK);
IK = p.gK.*nK.^4.*(V - p.vK);
IL = p.gL.*(V - p.vL);
IH = p.gH.*hH.*(V - p.vH);
INaP = p.gNaP.*mNaP_inf*hNaP*(V - p.vNa);
ILVA = p.gLVA.*mLVA.^2.*hLVA.*(V - p.vCa);
IHVA = p.gHVA*mHVA_inf*(V - p.vCa);
IBK = p.gBK*wBK*(V - p.vK);

%% get ORN input for t
Input = p.intrace(floor(t / p.ORNsamplingfactor)+1);


xdot = zeros(size(x));
xdot(1) = -(INa + IK + ILVA + IH + INaP + IL + IHVA + IBK + INew - Input)./p.C;
xdot(2) = (nK_inf-nK)/nK_tau;
xdot(3) = (hNaP_inf-hNaP)/hNaP_tau;
xdot(4) = (hH_inf-hH)./hH_tau;
xdot(5) = (mLVA_inf-mLVA)/mLVA_tau;
xdot(6) = (hLVA_inf-hLVA)/hLVA_tau;
xdot(7) = (wBK_inf - wBK)/wBK_tau;
xdot(8) = -p.Ca_buffer*10*(ILVA + IHVA)/(p.Ca_z*p.F*p.d) + (p.Ca0 - Ca)/p.tau_Ca;
xdot(9) = (nNew_inf - nNew)./nNew_tau;

end

function data = integrator(odefun, odeopts, ics, PARS, ORNtrace, ORNsamplingrate)

    ics = ics(:);

    PARS.intrace = ORNtrace;
    PARS.ORNsamplingfactor = 1000 / ORNsamplingrate;

    %% Integration tolerances
    atol = 1e-6;
    rtol = 1e-6;
    options = odeset('reltol',rtol,'abstol',atol,'InitialStep',0.5);
    options = odeset(options, odeopts); %override odeopts

    %% time span
    tstart = 0;
    tend = (length(ORNtrace) - 1) * PARS.ORNsamplingfactor;

    %% integrate
    tout = tstart;
    xout = ics';
    teout = [];
    ieout = [];

    tic
    while tstart < tend
        [t,x,te,xe,ie] = ode15s(odefun,[tstart tend],ics,options,PARS);

        nt = length(t);
        tout = [tout; t(2:nt)];
        xout = [xout; x(2:nt,:)];
        tstart = t(nt);

        if ~isempty(te)
            teout = [teout; te];
            ieout = [ieout; ie];

            ics = xe(end,:);

            options = odeset(options,'InitialStep',t(nt)-t(nt-1)); %use most recent timestep
        end
    end

    %%%%%%%%%%%%%calculate currents%%%%%%
    % 1 = transient sodium
    % 2 = fast potassium
    % 3 = leak
    % 4 = persistent sodium
    % 5 = hyperpolarization activated
    % 6 = LVA calcium
    % 7 = HVA calcium
    % 8 = large conductance potassium
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    current(:,1) = PARS.gNa.*(1-xout(:,2)).*(1./(1+exp((xout(:,1)-PARS.theta_mNa)./PARS.sigma_mNa))).^3.*(xout(:,1)-PARS.vNa);
    current(:,2) = PARS.gK.*xout(:,2).^4.*(xout(:,1)-PARS.vK);
    current(:,3) = PARS.gL.*(xout(:,1)-PARS.vL);
    current(:,4) = PARS.gNaP.*(1./(1+exp((xout(:,1)-PARS.theta_mNaP)./PARS.sigma_mNaP))).*xout(:,3).*(xout(:,1)-PARS.vNa);
    current(:,5) = PARS.gH.*xout(:,4).*(xout(:,1)-PARS.vH);
    current(:,6) = PARS.gLVA.*xout(:,5).^2.*xout(:,6).*(xout(:,1)-PARS.vCa);
    current(:,7) = PARS.gHVA.*1./(1+exp(-(xout(:,1) + 20)./9)).*(xout(:,1)-PARS.vCa);
    current(:,8) = PARS.gBK.*xout(:,7).*(xout(:,1)-PARS.vK);

    wBK_p = 2.9 + 6.3.*exp(-360.*xout(:,8));
    wBK_s = -25.3 + 107.5.*exp(-120.*xout(:,8));
    wBK_f = 1./(10.*(exp(-(xout(:,1) + 100 - wBK_s)./63.6)+exp(-150-(xout(:,1) + 100 - wBK_s))./63.6)) - 5.2;
    current(:,9) = -(wBK_p - 1).*(wBK_f - 0.2)./0.8 + PARS.wBK_base;

    data = struct('T', tout, 'X', xout, 'events', teout, 'which', ieout, 'current', current);
    toc

end


function [value,isterminal,direction] = spikedetect(~,y,~)
    value(1) = y(1); % spikes
    value(2) = y(2) - 0.25; % burst start
    value(3) = y(2) - 0.25; % burst end
    isterminal(1) = 0;
    isterminal(2) = 0;
    isterminal(3) = 0;
    direction(1) = 1;
    direction(2) = 1;
    direction(3) = -1;
end
