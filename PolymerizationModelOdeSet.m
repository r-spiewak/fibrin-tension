function dots = PolymerizationModelOdeSet(t,y,dotParams)
    %% Notes
    % f_A    = concentration of fibrinogen monomers
    % f      = concentration of fibrin monomers
    % f_n    = concentration of protofibrils
    % f_r    = concentration of fibrin fibers
    % f_nTot = total protofibrils in fibers
    % k_A    = rate constant: fibrinogen to fibrin monomer
    % k_pi   = rate constant: protofibril formation
    % k_pg   = rate constant: protofibril growth
    % k_fi   = rate constant: fiber formation
    % k_fg   = rate constant: fiber growth
    % m      = f_nTot/f_r, average number of protofibrils per fiber cross-section
    
    fLogE = dotParams{1};
    k_A  = dotParams{2};
    k_pi = dotParams{3};
    k_pg = dotParams{4};
    k_fi = dotParams{5};
    k_fg = dotParams{6};
    minLength = dotParams{7};
    k_fA = dotParams{8};
    if minLength<2
        minLength=2;
        warning('minLength must be >= 2. Setting minLength=2.');
    end
    
    f_A    = y(1);
    f = zeros(1,minLength);
    for ii=1:minLength
        f(ii)   = y(ii+1);
    end
    f_n    = y(minLength+2);
    f_r    = y(minLength+3);
    f_nTot = y(minLength+4);
    c_fn   = y(minLength+5);
    c_fr   = y(minLength+6);
    
    df_A    = -k_A*f_A;
    df = zeros(1,minLength);
    sumf = sum(f);
    for jj=1:minLength
        df(jj) = -f(jj).*f(jj) - f(jj).*sumf;
        for ii=1:floor(jj/2)
            df(jj) = df(jj) + f(ii).*f(jj-ii);
        end
        df(jj) = k_pi*df(jj) - k_pg*f_n.*f(jj);
    end
    df(1) = df(1) + k_A*f_A;
    df_n = 0;
    for jj=1:floor((minLength+1)/2)
        temp = 0;
        for ii=minLength+1-jj:minLength
            temp = temp + f(ii);
        end
        temp = temp.*(f(jj)+f(minLength+1-jj));
        df_n = df_n + temp;
    end
    df_n = k_pi*df_n - 2*k_fi*f_n.*f_n - k_fg*f_r.*f_n;
    df_r    = k_fi*f_n.*f_n - k_fA*f_r.*f_r;
    df_nTot = 2*k_fi*f_n.*f_n + k_fg*f_r.*f_n + k_fA*f_r.*f_r;
    dc_fn = 0;
    for ii=1:minLength
        temp = 0;
        for jj=ii:floor((minLength+ii)/2)
            temp = temp + f(jj).*f(minLength+ii-jj);
        end
        dc_fn = dc_fn +(minLength+ii)*temp;
    end
    dc_fn = k_pi*dc_fn + k_pg*f_n.*sumf - k_fi*f_n.*c_fn - k_fg*f_r.*c_fn;
    dc_fr   = 2*k_fi*f_n.*c_fn + k_fg*f_r.*c_fn + k_fA*f_r.*f_r;
    
    dots = zeros(minLength+6,1);
    dots(1) = df_A;
    for ii=2:minLength+1
        dots(ii) = df(ii-1);
    end
    dots(minLength+2) = df_n;
    dots(minLength+3) = df_r;
    dots(minLength+4) = df_nTot;
    dots(minLength+5) = dc_fn;
    dots(minLength+6) = dc_fr;
    
    
    if(toc>0.1)
        fprintf(fLogE, '\t%s: %s%6.3f\r\n', datestr(now,0),'ODE t=',t);
        tic;
    end
end