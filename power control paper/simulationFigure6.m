
close all; clear;

%% Simulation setup
nbrOfSetups       = 10;     % raise for final
nbrOfRealizations = 50;     % raise for final

L = 100;         % APs
N = 4;           % antennas/AP
K = 40;          % UEs
tau_c = 200; 
tau_p = 10;
preLogFactor = (tau_c - tau_p)/tau_c;

p_ul   = 100;    % mW, UL power for channel estimation / virtual UL
rho_tot = 200;   % mW, total DL power budget per AP

% Storage (only what we plot)
SE_DL_equal  = zeros(K, nbrOfSetups);
SE_DL_fpa    = zeros(K, nbrOfSetups);   % FPA: nu=-0.5, kappa=0.5
SE_DL_mmf    = zeros(K, nbrOfSetups);
SE_DL_sumSE  = zeros(K, nbrOfSetups);

%% Main loop
for n = 1:nbrOfSetups
    disp(['Setup ' num2str(n) ' / ' num2str(nbrOfSetups)]);
    
    % --- Generate single setup and slice
    [gainAll, Rall, pilotAll, Dall] = generateSetup(L, K, N, tau_p, 1);
    gainOverNoisedB = gainAll(:,:,1);         % L×K
    R               = Rall(:,:,:,:,1);        % N×N×L×K
    pilotIndex      = pilotAll(:,1);          % K×1
    D               = Dall(:,:,1);            % L×K

    % --- UL power vector (per UE)
    p_vec = p_ul * ones(K,1);

    % --- MMSE channel estimation
    [Hhat, H, ~, C] = functionChannelEstimates(R, nbrOfRealizations, L, K, N, tau_p, pilotIndex, p_vec);

    % --- Expectation terms for DL with P-MMSE (UL–DL duality)
    [signal_P_MMSE, signal2_P_MMSE, scaling_P_MMSE] = ...
        functionComputeExpectations(Hhat, H, D, C, nbrOfRealizations, N, K, L, p_vec);

    % ---- (7.13)–(7.15) terms
    bk = abs(diag(signal_P_MMSE)).^2;      % K×1
    ck = signal2_P_MMSE.';                 % K×K
    ck = ck - diag(bk);                    % remove desired

    % ---- Normalize by expected precoder norm (7.16)
    sigma2 = sum(scaling_P_MMSE, 1).';     % K×1
    sigma2 = max(sigma2, eps);
    bk = bk ./ sigma2;
    % row-wise divide: each row k of ck divided by sigma2(k)
    ck = ck ./ repmat(sigma2, 1, K);

    % ---- Portion scaling (per-AP power split of each UE’s precoder)
    denom = max(sum(scaling_P_MMSE, 1), eps);         % 1×K
    portionScaling = scaling_P_MMSE ./ denom;         % L×K

    % ---- DL power allocations for Fig. 7.2(b)
    rho_equal = (rho_tot / tau_p) * ones(K,1);
    % FPA with (nu=-0.5, kappa=0.5)
    rho_fpa = functionCentralizedPowerAllocation(K, gainOverNoisedB, D, rho_tot, portionScaling, -0.5, 0.5);

    % ---- Downlink SEs (Theorem 6.1)
    SE_DL_equal(:,n) = preLogFactor * log2(1 + bk .* rho_equal ./ (ck' * rho_equal + 1));
    SE_DL_fpa(:,n)   = preLogFactor * log2(1 + bk .* rho_fpa   ./ (ck' * rho_fpa   + 1));

    % ---- Max-min fair and Sum-SE allocations (centralized)
    SE_DL_mmf(:,n)   = functionDownlinkSE_maxmin(bk, ck, portionScaling, preLogFactor, K, rho_tot);
    SE_DL_sumSE(:,n) = functionDownlinkSE_sumSE (bk, ck, portionScaling, preLogFactor, L, K, rho_tot, tau_p);
end

%% ============================ Plot 7.2(b) ============================
figure; hold on; box on; set(gca,'fontsize',16);
plot(sort(SE_DL_equal(:)), linspace(0,1,K*nbrOfSetups), 'k-',  'LineWidth',2);
plot(sort(SE_DL_fpa(:)),   linspace(0,1,K*nbrOfSetups), 'k:',  'LineWidth',2);
plot(sort(SE_DL_mmf(:)),   linspace(0,1,K*nbrOfSetups), 'b-.', 'LineWidth',2);
plot(sort(SE_DL_sumSE(:)), linspace(0,1,K*nbrOfSetups), 'r--', 'LineWidth',2);
xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Equal','FPA','MMF','SumSE'}, 'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);
title('Figure 7.2(b) — DL P-MMSE');
