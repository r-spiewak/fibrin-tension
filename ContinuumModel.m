function [ls,le,Ft,T11,tanAlphaS,tanAlphaE,LL,lSt] = ContinuumModel(fLog,rm,a0,k0,tanAlpha0,K,Kb,tol,R,l,phiS,nuType,trackProgress)
    % nuType == 1 % Fixed number density (nd)
    % nuType == 2 % Fixed solid volume fraction in free fiber config
    % nuType == 3 % Fixed solid volume fraction in network config (vf)
    % 
    syms lambdaStar;
    
    lSt = 0.*R;
    ls  = 0.*R;
    le=0.*R;
    T11 = 0.*R;
    Ft = 0.*R;
    tanAlphaE=0.*R;
    tanAlphaS = sqrt(log(R./rm)./(k0.*(R-rm))-1);
    lsR = tanAlphaS./tanAlpha0;
    if lsR(1)==-Inf||lsR(1)==Inf||isnan(lsR(1))
        lsR(1)=0;
    end
    clear textprogressbar
    if(trackProgress==1)
        tstart=tic;
        tfin=length(R);
        textprogressbar('Continuum 8-Chain Model time progress: ');
    end
    for ii = 1:numel(R)
        if lsR(ii)==0
            continue
        end
        f = K.*(lambdaStar.^2./lsR(ii)-1);
        lambda_ce = sqrt((1+2*(lambdaStar.*lsR(ii)).^2)./(3.*lsR(ii).^2)); % same as lambda^e
        lambda_1e = 1./lsR(ii);
        L=l./(lambda_ce.*lsR(ii));
        if nuType==1 %fixed number density
            nu=3*sqrt(3)/L.^3;
        elseif nuType==2 %fixed solid volume fraction in free fiber config
            nu=phiS./(pi*R(ii).^2.*lsR(ii).*L);
        elseif nuType==3 %fixed solid volume fraction in network config
            nu=phiS./(pi*R(ii).^2.*lsR(ii).*lambda_ce.*L);
        else
            error('Invalid nuType %s.',nuType);
        end
        tanAlpha = lambda_ce.*tanAlphaS(ii);
        hyp = sqrt(1+tanAlpha.^2);
        sinAlpha = tanAlpha./hyp;
        cosAlpha = 1/hyp;
        F = (2*pi/a0)*Kb*sinAlpha.*(-log(R(ii)./rm).*(cosAlpha).^2+k0.*(R(ii)-rm));
        C1 = nu*L/3*k0*(R(ii)-rm).*2.*pi.*Kb./a0.*tanAlphaS(ii);
        C2 = nu*L/3*log(R(ii)./rm).*2.*pi.*Kb./a0.*tanAlphaS(ii);
        C3 = lsR(ii)*K;
        C4 = subs(1+1/(3*lsR(ii).^2).*((log(R(ii)./rm))./(k0.*(R(ii)-rm))-1));
        C5 = 2/3*((log(R(ii)./rm))./(k0*(R(ii)-rm))-1);
        A0 = C3.^2.*C4.^3-(C1.*C4-C2).^2;
        A1 = 3*C3.^2.*C4.^2.*C5-2*C3.*C4.^3.*K-2*(C1.*C4-C2).*C1.*C5;
        A2 = 3*C3.^2.*C4.*C5.^2+C4.^3.*K.^2-6*C3.*C4.^2.*C5.*K-(C1.*C5).^2;
        A3 = C3.^2.*C5.^3+3*C4.^2.*C5.*K.^2-6*C3.*C4.*C5.^2.*K;
        A4 = 3*C4.*C5.^2.*K.^2-2*C3.*C5.^3.*K;
        A5 = C5.^3.*K.^2;
        tenthDegreePolynomial = A5*lambdaStar.^10 + A4*lambdaStar.^8 + A3*lambdaStar.^6 + A2*lambdaStar.^4 + A1*lambdaStar.^2 + A0 == 0;
        S10=vpasolve(tenthDegreePolynomial);
        clear S
        S = find(S10>0&imag(S10)==0);
        if isempty(S)
            message = sprintf('No viable solutions found for polynomial on iteration %d.',ii);
            warning(message)
            fprintf(fLog, '%s: Warning: %s.\r\n', datestr(now,0),message);
            %dbstop
            %keyboard
            continue
        else
            res = zeros(size(S));
            accept = zeros(size(S));
            for jj=1:numel(S)
                res(jj) = double(subs(nu*L ./(3*lambda_ce).*F+lsR(ii).*f,[lambdaStar],[S10(S(jj))]));
                if abs(res(jj))<tol
                    accept(jj) = 1;
                end
            end
            if sum(accept)==1
                lSt(ii) = S10(sum(S.*accept));
            elseif sum(accept)==0
                message = sprintf('No viable solutions that make T22=0 found for polynomial on iteration %d.',ii);
                warning(message)
                fprintf(fLog, '%s: Warning: %s.\r\n', datestr(now,0),message);
                lSt(ii) = 0;
            else %multiple acceptable solutions
                lSt(ii) = min(S10(setdiff(S.*accept,0)));
                message = sprintf('Multiple (%d) viable solutions that make T22=0 found for polynomial on iteration %d. Using the smallest solution: %f.',sum(accept),ii,lSt(ii));
                warning(message)
                fprintf(fLog, '%s: Warning: %s.\r\n', datestr(now,0),message);
            end
        end
        ls(ii)=lsR(ii);
        le(ii)=subs(lambda_ce,lambdaStar,lSt(ii));
        tanAlphaE(ii) = double(subs(tanAlpha,[lambdaStar],[lSt(ii)]));
        % Force on single fiber
        Ft(ii) = subs(F,[lambdaStar],[lSt(ii)]);
        % Piola Stress
        T11(ii) = subs(nu*L*lambda_1e./(3.*lambda_ce).*F +lsR(ii).^2.*lambdaStar.^2.*f,[lambdaStar],[lSt(ii)]);
        if (trackProgress == 1)
            progress=100*ii/tfin;
            textprogressbar(progress);
        end
    end
    if (trackProgress==1)
        tfin=[];
        textprogressbar('');
        display([ '   Continuum 8-Chain Model time completion time: ' num2str(toc(tstart))]);
        tstart=[];
    end
    LL=l./(le(end).*ls(end));
end