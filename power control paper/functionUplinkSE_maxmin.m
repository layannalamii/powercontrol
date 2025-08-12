function SE = functionUplinkSE_maxmin(bk,ck,sigma2,preLogFactor,K,pmax)

%Implement Algorithm 7.1

%Initialize iteration counter to zero
iter = 0;
%Initialize the power control coefficients to full power
eta = pmax*ones(K,1);
%Compute the denominator in (7.1) for all UEs
denominator = ck'*eta+sigma2;
%Compute SINRs in (7.1) for all UEs
SINR = eta.*bk./denominator;

while max(SINR)-min(SINR)>0.01 %the condition in Line 2 of Algorithm 7.1 with solution accuracy 0.01
    %Increase iteration counter by one
    iter = iter+1;
    
    eta = denominator./bk; %Line 3 of Algorithm 7.1
    eta = eta*pmax/max(eta); %Line 4 of Algorithm 7.1
    
    %Update the denominator in (7.1) for all UEs
    denominator = ck'*eta+sigma2;
    %Update SINRs in (7.1) for all UEs
    SINR = eta.*bk./denominator;
   
end
%Compute SEs
SE = preLogFactor*log2(1+SINR);

