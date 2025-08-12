% Complexity Analysis for Scalable Cell-Free Massive MIMO (Log-Scale, Selected Detectors)
close all; clear;

%% Simulation parameters
nbrOfSetups = 2;          % Increase for publication quality
L = 100;
N = 4;
Kvec = 20:20:100;         % User counts
tau_c = 200;
tau_p = 10;
maxK = max(Kvec);

% Preallocate
D_tot = zeros(L, maxK, length(Kvec), nbrOfSetups);

%% Generate D matrices
disp('Generating DCC matrices ...');
for n = 1:nbrOfSetups
    disp(['Setup ' num2str(n) ' of ' num2str(nbrOfSetups)]);
    for nK = 1:length(Kvec)
        K = Kvec(nK);
         [gainOverNoisedB,R,pilotIndex,D,APpositions,UEpositions,distances] = generateSetup(L,K,N,tau_p,1);
        D_tot(:, 1:K, nK, n) = D;
    end
end

%% Complexity computation
MMSE_all = zeros(1, length(Kvec));
MMSE_DCC = zeros(1, length(Kvec));
PMMSE_DCC = zeros(1, length(Kvec));
LMMSE_all = zeros(1, length(Kvec));
LMMSE_DCC = zeros(1, length(Kvec));
LPMMSE_DCC = zeros(1, length(Kvec));

for nK = 1:length(Kvec)
    K = Kvec(nK);

    % MMSE (ALL) - Table 5.1
    MMSE_all(nK) = (N*tau_p+N^2)*K*L + ((N*L)^2+N*L)/2*K + (N*L)^2*K + ((N*L)^3-N*L)/3;

    % L-MMSE (ALL) - Table 5.3
    LMMSE_all(nK) = (N*tau_p+N^2)*K*L + (N^2+N)/2*K*L + N^2*K*L + (N^3-N)/3*L;

    % Now, for DCC/user-centric methods:
    MMSE_DCC_sum = 0;
    PMMSE_DCC_sum = 0;
    LMMSE_DCC_sum = 0;
    LPMMSE_DCC_sum = 0;

    for n = 1:nbrOfSetups
        Dn = D_tot(:, 1:K, nK, n);

        % Number of unique APs used in this setup
        L_used = sum(sum(Dn,2) >= 1);

        % MMSE (DCC) - FULL formula (Table 5.1, replace L with L_used)
        MMSE_DCC_sum = MMSE_DCC_sum + ...
            (N*tau_p+N^2)*K*L_used + ...
            ((N*L_used)^2+N*L_used)/2*K + ...
            (N*L_used)^2*K + ((N*L_used)^3-N*L_used)/3;

        % L-MMSE (DCC) - Table 5.3, replace L with L_used
        LMMSE_DCC_sum = LMMSE_DCC_sum + ...
            (N*tau_p+N^2)*K*L_used + (N^2+N)/2*K*L_used + N^2*K*L_used + (N^3-N)/3*L_used;

        % For user-centric (P-MMSE/LP-MMSE), sum over users
        for k = 1:K
            servingAPs = find(Dn(:, k)==1);
            La = length(servingAPs);
            servedUEs = find(sum(Dn(servingAPs, :),1) >= 1);
            Bk = length(find(sum(Dn(:, servedUEs),2) >= 1));

            % P-MMSE (DCC)
            PMMSE_DCC_sum = PMMSE_DCC_sum + ...
                (N*tau_p+N^2)*Bk + ((N*Bk)^2+N*Bk)/2 + (N*La)^2 + ((N*La)^3-N*La)/3;

            % LP-MMSE (DCC)
            LPMMSE_DCC_sum = LPMMSE_DCC_sum + ...
                (N*tau_p+N^2)*La + (N^2+N)/2*La + N^2*La + (N^3-N)/3*La;
        end
    end

    MMSE_DCC(nK) = MMSE_DCC_sum / nbrOfSetups;
    PMMSE_DCC(nK) = PMMSE_DCC_sum / nbrOfSetups;
    LMMSE_DCC(nK) = LMMSE_DCC_sum / nbrOfSetups;
    LPMMSE_DCC(nK) = LPMMSE_DCC_sum / nbrOfSetups;
end

%% Plot: Centralized Detectors (log-scale)
figure; hold on; box on; grid on;
set(gca, 'fontsize', 16, 'YScale', 'log');
semilogy(Kvec, MMSE_all/L, 'kd-', 'LineWidth', 2, 'MarkerFaceColor','k');
semilogy(Kvec, MMSE_DCC/L, 'bs--', 'LineWidth', 2, 'MarkerFaceColor','b');
semilogy(Kvec, PMMSE_DCC/L, 'ro-.', 'LineWidth', 2, 'MarkerFaceColor','r');
xlabel('Number of UEs, $K$', 'Interpreter', 'Latex');
ylabel('Number of complex multiplications per AP', 'Interpreter', 'Latex');
legend({'MMSE (All)','MMSE (DCC)','P-MMSE (DCC)'}, ...
    'Interpreter', 'Latex', 'Location', 'SouthWest');
title('\bf Centralized operation', 'Interpreter', 'Latex');
set(gca, 'XTick', Kvec);

%% Plot: Distributed Detectors (log-scale)
figure; hold on; box on; grid on;
set(gca, 'fontsize', 16, 'YScale', 'log');
semilogy(Kvec, LMMSE_all/L, 'kd-', 'LineWidth', 2, 'MarkerFaceColor','k');
semilogy(Kvec, LMMSE_DCC/L, 'bs--', 'LineWidth', 2, 'MarkerFaceColor','b');
semilogy(Kvec, LPMMSE_DCC/L, 'mo-.', 'LineWidth', 2, 'MarkerFaceColor','m');
xlabel('Number of UEs, $K$', 'Interpreter', 'Latex');
ylabel('Number of complex multiplications per AP', 'Interpreter', 'Latex');
legend({'L-MMSE (All)','L-MMSE (DCC)','LP-MMSE (DCC)'}, ...
    'Interpreter', 'Latex', 'Location', 'SouthWest');
title('\bf Distributed operation', 'Interpreter', 'Latex');
set(gca, 'XTick', Kvec);
