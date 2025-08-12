function [signal_P_MMSE, signal2_P_MMSE, scaling_P_MMSE] = ...
    functionComputeExpectations(Hhat,H,D,C,nbrOfRealizations,N,K,L,p)

eyeN = eye(N);
PowMat = diag(p);

% Scale C by user powers (uplink)
Cp = zeros(size(C));
for k=1:K
    Cp(:,:,:,k) = p(k)*C(:,:,:,k);
end

signal_P_MMSE  = zeros(K,K);
signal2_P_MMSE = zeros(K,K);
scaling_P_MMSE = zeros(L,K);

for n=1:nbrOfRealizations
    
    for k = 1:K
        servingAPs = find(D(:,k)==1);
        La = length(servingAPs);
        servedUEs = sum(D(servingAPs,:),1)>=1;

        Hallj_active    = zeros(N*La,K);
        Hhatallj_active = zeros(N*La,K);
        Cp_blk_partial  = zeros(N*La,N*La);

        for l = 1:La
            idx = servingAPs(l);
            Hallj_active((l-1)*N+1:l*N,:)    = reshape(H((idx-1)*N+1:idx*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((idx-1)*N+1:idx*N,n,:),[N K]);
            Cp_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(Cp(:,:,idx,servedUEs),4);
        end

        Hphat = Hhatallj_active * sqrt(PowMat);

        % P‑MMSE combiner
        w = ((Hphat(:,servedUEs)*Hphat(:,servedUEs)' + Cp_blk_partial + eye(La*N)) \ Hphat(:,k)) * sqrt(p(k));

        tempor = Hallj_active' * w;
        signal2_P_MMSE(:,k) = signal2_P_MMSE(:,k) + abs(tempor).^2 / nbrOfRealizations;
        signal_P_MMSE(:,k)  = signal_P_MMSE(:,k)  + tempor          / nbrOfRealizations;

        for l=1:La
            w_l = w((l-1)*N+1:l*N,:);
            scaling_P_MMSE(servingAPs(l),k) = scaling_P_MMSE(servingAPs(l),k) + sum(abs(w_l).^2) / nbrOfRealizations;
        end
    end
end
end
