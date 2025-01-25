function [ETscm,ETgafr,ETcicm] = ExecTimeMaxAlgTridFunc(J,epsi,Nr)
%
% This function computes the average time to reach an epsi-optimal solution of 
% a trid function when three SCM, GAFR and CICM methods are respectively used
%
% -------------------------------- INPUT -------------------------------- % 
% J: number of variables (i.e., dimension of the trid function)
% epsi: acceptable error to the maximum
% Nr: the number of repetitions to compute the average quantities
%
% ------------------------------- OUTPUT -------------------------------- %
% ETscm, ETgafr, ETcicm: average execution time of SCM, GAFR and CICM methods
%
% ------------------------------ EXAMPLE -------------------------------- %
% J = 100; epsi = 1e-6; Nr = 10; 
% [ETscm,ETgafr,ETcicm] = ExecTimeMaxAlgTridFunc(J,epsi,Nr)
%
% ---------------------------- MAIN CODE -------------------------------- %

% compution of the theoretical optimal value of the trid function
Fopt = J*(J+4)*(J-1)./6;

% random choice of initial vector and compution of associated function value
% pre-allocation
vFint = zeros(1,Nr); vX0 = zeros(J,Nr);
% loop
for r = 1:Nr
    % choose randomly initial vector over [-J^2,J^2]^J
    vX0(:,r) = -J^2 + 2*J^2*rand(J,1);
    % compute initial value vector of the trid function
    vFint(r) = -(vX0(:,r)-1)'*(vX0(:,r)-1) + (vX0(2:end,r))'*vX0(1:(end-1),r);
end

% Pre-allocation for vector of execute times, optimal value and optima
vTscm = zeros(1,Nr); vTgafr = vTscm; vTcicm = vTscm;

% ------------------------- Start of SCM method ------------------------- %

% Main loop
for r = 1:Nr
        
    % start a stopwatch timer
    tstartscm = cputime;
   
    % compute the initial residual of function value
    rscm = Fopt - vFint(r);
    
	% set initial value of Xoptscm
	Xoptscm = vX0(:,r);
    
    % repeat the loop when the condition rscm > epsi is till true
    while rscm > epsi
               
        % pre-allocation for vector of a_j^(k) 
        vaj = zeros(J,1);
        
        % computation of a_j^(k) following condition of j
        vaj(1) = Xoptscm(2)/2 - Xoptscm(1) + 1;
        vaj(2:(J-1)) = ( Xoptscm((2:(J-1))-1)+Xoptscm((2:(J-1))+1) )/2 - Xoptscm((2:(J-1))) + 1;
        vaj(J) =  Xoptscm(J-1)/2 - Xoptscm(J) + 1;
               
        % computation of a_0^(k)
        a0 = ( -2*(vaj'*(Xoptscm - 1)) + (vaj(1:(end-1)))'*Xoptscm(2:end) + (vaj(2:end))'*Xoptscm(1:(end-1)) ) /...
               ( 2*(vaj'*vaj) - 2*(vaj(1:(end-1)))'*vaj(2:end) );
        
        % update value of Xoptscm
        Xoptscm = Xoptscm + a0*vaj;
        
        % update the function value
        Foptscm = - (Xoptscm-1)'*(Xoptscm-1) + (Xoptscm(2:end))'*Xoptscm(1:(end-1));
        
        % update the residual of function value
        rscm = Fopt - Foptscm; 
        
    end
    
    % Read elapsed time from stopwatch
    vTscm(r) = cputime - tstartscm;

end

% compute the average execution time
ETscm = mean(vTscm);

% -------------------------- end of SCM method -------------------------- %

% ------------------------ Start of GAFR method ------------------------- %

% Main loop
for r = 1:Nr
        
    % start a stopwatch timer
    tstartgafr = cputime;
   
    % compute the initial residual of function value
    rgafr = Fopt - vFint(r);
    
	% set initial value of Xoptgafr
	Xoptgafr = vX0(:,r);
    
    % repeat the loop when the condition rgafr > epsi is till true
    while rgafr > epsi
               
        % pre-allocation for vector of gradient 
        vgrad = zeros(J,1);
        
        % computation of gradient following condition of j
        vgrad(1) = -2*(Xoptgafr(1) - 1) + Xoptgafr(2);
        vgrad(2:(J-1)) = -2*(Xoptgafr(2:(J-1)) - 1) + Xoptgafr((2:(J-1))-1) + Xoptgafr((2:(J-1))+1);
        vgrad(J) = -2*(Xoptgafr(J) - 1) + Xoptgafr(J-1);
               
        % computation of a_0^(k)
        a0 = ( -2*(vgrad'*(Xoptgafr - 1)) + (vgrad(1:(end-1)))'*Xoptgafr(2:end) + (vgrad(2:end))'*Xoptgafr(1:(end-1)) ) /...
               ( 2*(vgrad'*vgrad) - 2*(vgrad(1:(end-1)))'*vgrad(2:end) );
        
        % update value of Xoptgafr
        Xoptgafr = Xoptgafr + a0*vgrad;
        
        % update the function value         
        Foptgafr = - (Xoptgafr-1)'*(Xoptgafr-1) + (Xoptgafr(2:end))'*Xoptgafr(1:(end-1));
        
        % update the residual of function value
        rgafr = Fopt - Foptgafr; 
        
    end
    
    % Read elapsed time from stopwatch
    vTgafr(r) = cputime - tstartgafr;

end

% compute the average execution time
ETgafr = mean(vTgafr);

% -------------------------- end of GAFR method ------------------------- %

% ------------------------ Start of CICM method ------------------------- %

% Main loop
for r = 1:Nr
        
    % start a stopwatch timer
    tstartcicm = cputime;
   
    % compute the initial residual of function value
    rcicm = Fopt - vFint(r);
    
	% set initial value of Xoptcicm
	Xoptcicm = vX0(:,r);
    
    % repeat the loop when the condition rcicm > epsi is till true
    while rcicm > epsi
        
        % update value of Xoptcicm loop
        for j = 1:J
            if j == 1
                Xoptcicm(j) = Xoptcicm(j+1)/2 + 1;
            elseif j == J
                Xoptcicm(j) = Xoptcicm(j-1)/2 + 1;
            else
                Xoptcicm(j) = ( Xoptcicm(j-1) + Xoptcicm(j+1) )/2 + 1;
            end
        end
                       
        % update the function value
        Foptcicm = - (Xoptcicm-1)'*(Xoptcicm-1) + (Xoptcicm(2:end))'*Xoptcicm(1:(end-1));
        
        % update the residual of function value
        rcicm = Fopt - Foptcicm; 
        
    end
    
    % Read elapsed time from stopwatch
    vTcicm(r) = cputime - tstartcicm;

end

% compute the average execution time
ETcicm = mean(vTcicm);

% -------------------------- end of CICM method ------------------------- %

end
