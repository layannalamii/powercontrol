function [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p)

%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
H = (randn(L*N,nbrOfRealizations,K)+1i*randn(L*N,nbrOfRealizations,K));

%Go through all channels and apply the spatial correlation matrices
for l = 1:L
    for k = 1:K
        Rsqrt = sqrtm(R(:,:,l,k));
        H((l-1)*N+1:l*N,:,k) = sqrt(0.5)*Rsqrt*H((l-1)*N+1:l*N,:,k);
    end
end

%% Perform channel estimation

eyeN = eye(N);
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,L,tau_p) + 1i*randn(N,nbrOfRealizations,L,tau_p));
Hhat = zeros(L*N,nbrOfRealizations,K);

if nargout>2
    B = zeros(size(R));
end
if nargout>3
    C = zeros(size(R));
end

%Go through all APs
for l = 1:L
    %Go through all pilots
    for t = 1:tau_p
        usersWithPilot = find(t==pilotIndex); % UEs using pilot t
        % --- Form processed pilot signal (sum for all UEs using this pilot, weighted by sqrt of power) ---
        yp = zeros(N, nbrOfRealizations);
        for idx = 1:numel(usersWithPilot)
            k_idx = usersWithPilot(idx);
            yp = yp + sqrt(p(k_idx))*tau_p*H((l-1)*N+1:l*N,:,k_idx);
        end
        yp = yp + sqrt(tau_p)*Np(:,:,l,t);

        % --- Form PsiInv: sum power-weighted R for all UEs using this pilot ---
        PsiInv = eyeN;
        for idx = 1:numel(usersWithPilot)
            k_idx = usersWithPilot(idx);
            PsiInv = PsiInv + p(k_idx)*tau_p*R(:,:,l,k_idx);
        end

        % --- Now estimate each user that uses this pilot ---
        for k = usersWithPilot'
            RPsi = R(:,:,l,k) / PsiInv;
            Hhat((l-1)*N+1:l*N,:,k) = sqrt(p(k))*RPsi*yp;
            if nargout>2
                B(:,:,l,k) = p(k)*tau_p*RPsi*R(:,:,l,k);
            end
            if nargout>3
                C(:,:,l,k) = R(:,:,l,k) - B(:,:,l,k);
            end
        end
    end
end
