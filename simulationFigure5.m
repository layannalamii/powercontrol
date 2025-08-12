% Empty workspace and close figures
close all; clear;

%% Define simulation setup
nbrOfSetups = 20;
nbrOfRealizations = 50;

L = 100; N = 4; K = 40;

tau_c = 200; tau_p = 10;
preLogFactor = (tau_c - tau_p)/tau_c;

% --- Power (we'll build per-UE vectors below)
p_common = 100;  % mW

% Prepare to save centralized P-MMSE results
SE_UL_PMMSE_full        = zeros(K,nbrOfSetups);
SE_UL_PMMSE_fractional  = zeros(K,nbrOfSetups);
SE_UL_PMMSE_fractional2 = zeros(K,nbrOfSetups);
SE_UL_PMMSE_maxmin      = zeros(K,nbrOfSetups);
SE_UL_PMMSE_sumSE       = zeros(K,nbrOfSetups);

for n = 1:nbrOfSetups
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    % -------------------------
    % CHANGED: use generateSetup #2
    % Option A: get all setups at once, then slice nth
    [gainAll, Rall, pilotAll, Dall] = generateSetup(L,K,N,tau_p,nbrOfSetups);
    gainOverNoisedB_n = gainAll(:,:,n);             % L×K
    R_n               = Rall(:,:,:,:,n);            % N×N×L×K
    pilotIndex_n      = pilotAll(:,n);              % K×1
    D_n               = Dall(:,:,n);                % L×K
    % (If you prefer, call generateSetup(L,K,N,tau_p,1) inside the loop and drop the “_n” indexing.)
    % -------------------------

    % CHANGED: estimator #2 expects per-UE powers
    % Build the three power-control vectors
    p_full        = p_common * ones(K,1);
    p_fractional  = zeros(K,1);
    p_fractional2 = zeros(K,1);

    % Large-scale gains (linear) for FPC
    gainOverNoise_lin = db2pow(gainOverNoisedB_n);

    % Compute FPC coefficients within each UE’s cooperation cluster (same logic as before)
    for k = 1:K
        servingAPs  = find(D_n(:,k)==1);
        servedUEs   = find(sum(D_n(servingAPs,:),1) >= 1);
        norm_den    = 0;
        norm_den2   = 0;
        for i = servedUEs
            servingAPsi = find(D_n(:,i)==1);
            norm_den  = max(norm_den,  sqrt(sum(gainOverNoise_lin(servingAPsi,i))));
            norm_den2 = max(norm_den2, 1/sqrt(sum(gainOverNoise_lin(servingAPsi,i))));
        end
        p_fractional(k)  = p_common * sqrt(sum(gainOverNoise_lin(servingAPs,k))) / norm_den;
        p_fractional2(k) = p_common / sqrt(sum(gainOverNoise_lin(servingAPs,k))) / norm_den2;
    end

    % -------------------------
    % CHANGED: pass sliced R_n and per-UE power vectors to estimator #2
    [Hhat_full, H_full, ~, C_full] = functionChannelEstimates(R_n, nbrOfRealizations, L, K, N, tau_p, pilotIndex_n, p_full);
    [Hhat_fpc,  H_fpc,  ~, C_fpc ] = functionChannelEstimates(R_n, nbrOfRealizations, L, K, N, tau_p, pilotIndex_n, p_fractional);
    [Hhat_fpc2, H_fpc2, ~, C_fpc2] = functionChannelEstimates(R_n, nbrOfRealizations, L, K, N, tau_p, pilotIndex_n, p_fractional2);
    % -------------------------

    % Expectations for P-MMSE (unchanged API)
    [signal_full,  signal2_full,  scaling_full ] = functionComputeExpectations(Hhat_full, H_full, D_n, C_full,  nbrOfRealizations, N, K, L, p_full);
    [signal_fpc,   signal2_fpc,   scaling_fpc  ] = functionComputeExpectations(Hhat_fpc,  H_fpc,  D_n, C_fpc,   nbrOfRealizations, N, K, L, p_fractional);
    [signal_fpc2,  signal2_fpc2,  scaling_fpc2 ] = functionComputeExpectations(Hhat_fpc2, H_fpc2, D_n, C_fpc2,  nbrOfRealizations, N, K, L, p_fractional2);

    % Centralized P-MMSE SE calc
    bk_full  = abs(diag(signal_full)).^2;
    bk_fpc   = abs(diag(signal_fpc)).^2;
    bk_fpc2  = abs(diag(signal_fpc2)).^2;

    ck_full  = signal2_full  - diag(bk_full);
    ck_fpc   = signal2_fpc   - diag(bk_fpc);
    ck_fpc2  = signal2_fpc2  - diag(bk_fpc2);


    sigma2_full  = sum(scaling_full, 1).';
    sigma2_fpc   = sum(scaling_fpc,  1).';
    sigma2_fpc2  = sum(scaling_fpc2, 1).';

    SE_UL_PMMSE_full(:,n)        = preLogFactor * log2(1 + bk_full .*  p_full        ./ (ck_full'  * p_full        + sigma2_full ));
    SE_UL_PMMSE_fractional(:,n)  = preLogFactor * log2(1 + bk_fpc  .*  p_fractional  ./ (ck_fpc'   * p_fractional  + sigma2_fpc  ));
    SE_UL_PMMSE_fractional2(:,n) = preLogFactor * log2(1 + bk_fpc2 .*  p_fractional2 ./ (ck_fpc2'  * p_fractional2 + sigma2_fpc2));

    % Power-control optimizations on centralized model (unchanged)
    SE_UL_PMMSE_maxmin(:,n) = functionUplinkSE_maxmin(bk_full, ck_full, sigma2_full, preLogFactor, K, p_common);
    SE_UL_PMMSE_sumSE(:,n)  = functionUplinkSE_sumSE (bk_full, ck_full, sigma2_full, preLogFactor, K, p_common);
end

%% Plot P-MMSE CDF (unchanged)
figure; hold on; box on; set(gca,'fontsize',16);
plot(sort(SE_UL_PMMSE_full(:)),       linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_UL_PMMSE_fractional(:)), linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
ppp = plot(sort(SE_UL_PMMSE_fractional2(:)), linspace(0,1,K*nbrOfSetups),'k:o','LineWidth',2);
ppp.MarkerSize = 6; ppp.MarkerIndices = 1:ceil(K*nbrOfSetups/7):K*nbrOfSetups;
plot(sort(SE_UL_PMMSE_maxmin(:)),     linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);
plot(sort(SE_UL_PMMSE_sumSE(:)),      linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Full','FPC, $\upsilon=0.5$','FPC, $\upsilon=-0.5$','MMF','SumSE' },'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);
