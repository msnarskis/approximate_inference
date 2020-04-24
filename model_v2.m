function [resp] = model_v2(stim, cond, sig_t, sig_n, sig_v, sig_sa, sig_sv, pr_R, pr_C, noise_a, noise_v, nsamp)
% Return a matrix (n: locations, w: experimental W's) with probablity of
% answering R=1 (all on same side).
% match:  cond = 1 
% center: cond = 0
% NOTE: R determined from W --> all(w==1) || all(w==-1)

    trials = 100; % number of trial
    
    % size(stim) = (k: frames, n: location, W: corresponds to correct R)
    k = size(stim,1);
    n = size(stim,2);
    
    % Generate W, R vectors
    
    W = ones(2^k,k);
    for w = 1:2^k
        c = dec2bin(w-1,k);
        for i=1:k
            W(w,i) = W(w,i) * sign(str2num(c(i))-0.5);
        end
    end
    
    R = [1 zeros(1,2^k-2) 1];
    
    % init vars
    resp = zeros(n, size(stim,3));
    zero = zeros(k,n);
    
    % loop over data (:,:,w)
    for w = 1:size(stim,3)
        
        % init dat vars
        R_i = ones(n,1); % P(R=1)
        
        % loop over trials
        for t = 1:trials
            % noisify stimulus
            X_t =    (stim(:,:,w)+randn(k,n)*noise_a);
            X_n = -1*(stim(:,:,w)+randn(k,n)*noise_a);
            X_R =    (abs(stim(:,:,w)*cond+randn(k,n)*noise_v));
            X_L = -1*(abs(stim(:,:,w)*cond+randn(k,n)*noise_v));

            % calculate values of analytical variables
            a_tn = sig_n/(sig_t + sig_n);
            X_tn = X_t.*a_tn - X_n.*(1-a_tn);
            sig_tn = sig_t * a_tn;
            
            a_stn = sig_sa/(sig_tn + sig_sa);
            X_stn = X_tn*a_stn;
            sig_stn = sig_tn * a_stn;
            
            a_RL = 1/2;
            X_RL = X_R.*a_RL - X_L.*(1-a_RL);
            sig_RL = sig_v * a_RL;
            
            a_sRL = sig_sv/(sig_RL + sig_sv);
            X_sRL = X_RL*a_sRL;
            sig_sRL = sig_RL * a_sRL;
                        
            % calculations for R
            p_R0 = zeros(1,n);
            p_R1 = zeros(1,n);
            for wi = 1:size(W,1) % marginalize over W
                Wi = repmat(W(wi,:)',1,n);
                
                % W_i dependant variables
                a_tnRL = sig_sRL/(sig_stn + sig_sRL);
                X_tnRL = X_stn.*a_tnRL + Wi.*X_sRL.*(1-a_tnRL);
                sig_tnRL = sig_stn * a_tnRL;
                
                if R(wi)
                % P(X|R=1, W_i)
                    p_R1 = p_R1 + prod(...
                        (1-pr_C) * normcdf(zero, -Wi.*X_stn, sig_stn)...
                        .*normcdf(zero, -X_sRL, sig_sRL)...
                        + pr_C * normpdf(X_stn, Wi.*X_sRL, sig_stn + sig_sRL).*...
                        normcdf(zero, -Wi.*X_tnRL, sig_tnRL)...
                        ,1)/2;
                
                else
                    % P(X|R=0, W_i)
                    p_R0 = p_R0 + prod(...
                        (1-pr_C) * normcdf(zero, -Wi.*X_stn, sig_stn)...
                        .*normcdf(zero, -X_sRL, sig_sRL)...
                        + pr_C * normpdf(X_stn, Wi.*X_sRL, sig_stn + sig_sRL).*...
                        normcdf(zero, -Wi.*X_tnRL, sig_tnRL)...
                        ,1)/(2^k-2);
                end
            end
            
            % consolidate data
            p_R0 = p_R0 * (1-pr_R);
            p_R1 = p_R1 * pr_R;

            debug_p_r = p_R1 ./ (p_R0 + p_R1);
            
            p_Ri = sampling(p_R1 ./ (p_R0 + p_R1), nsamp);
            R_i = (R_i*(t-1) + p_Ri')./t;

        end % trials
        
        % consolidate data
        resp(:,w) = R_i
        
    end % data    
end


function [resp] = sampling(inp, nsamp)
    
    inp(isnan(inp)) = .5;
    
    if any(isnan(inp))
        ["num error:", 'w',w,'t',t]
    end
    
    % nsamples < Inf
    if nsamp~=Inf
        resp = binocdf(floor(nsamp/2), nsamp, inp, 'upper')...
        + .5*binopdf(nsamp/2, nsamp, inp);
    end

    % nsamples = Inf
    if nsamp==Inf
        resp = (inp > 0.5) + 0.5*(inp == 0.5);
    end
end
