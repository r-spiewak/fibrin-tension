function wkspacename = PiolaStressFromPolymerization(...
    savePic,haveExperimentalData,doPDESegmentLengthEquilibrium,trackProgress,...
    tolContinuum,fs,lw,linestyles,...
    xlimits,ylimits,...
    tstep,tlim,f_A0,f_n0,f_r0,f_nTot0,c_fn0,c_fr0,f0,...
    minLength,k_A,k_pi,k_pg,k_fi,k_fg,k_fA,...
    K,l,nu,rm,a0,Al,Kb,k0,tanAlpha0,mu,...
    nuType,phiS,...
    tMin,tFinal,tRestartSize,...
    Dt,lenDiffEnd,...
    scriptfignum,...
    restart,maxRestarts,restartWkspaceName,MainWkspaceName,maxArraySize,...
    filename...
    )
dbstop if error

tTotalStart=tic;
colorVector;

fLog = fopen(filename+"_logFile.txt",'w+');
if fLog == -1
    error('Cannot open log file.');
end

%% ODE System Parameters and Solve

fprintf(fLog, '%s: Loading input parameters for ODE system.\r\n', datestr(now,0));
savefiles = savePic; %1 to save the plots automatically, 0 to not save
fignum=1;
savefname = char(filename+"__Figure_"+fignum);

parameters = {
    savefname    %1
    savefiles
    tstep
    tlim
    f_A0         %5
    f0          
    f_n0
    f_r0
    f_nTot0
    k_A          %10
    k_pi
    k_pg
    k_fi
    k_fg
    xlimits      %15
    ylimits
    c_fn0
    c_fr0
    minLength
    k_fA         %20
    };
fprintf(fLog, '%s: Solving ODE system.\r\n', datestr(now,0));
[t,m,n,lCM,f_A,f,f_n,f_r,f_nTot,c_fn,c_fr] = PolymerizationModel(parameters);

fignum=fignum+1;
figname = filename+"__Figure_"+fignum; 
figure('Name',figname,'NumberTitle','on')
hold on
grid;
xlabel('Time (s)','FontSize',fs)
ylabel('fibrin per protofibril','FontSize',fs)
plot(t,n,'-','Color',colorVec(iiBlack,:),'linewidth',lw)
if(savePic==1)
    saveFig(figname);
end

fignum=fignum+1;
figname = filename+"__Figure_"+fignum; 
figure('Name',figname,'NumberTitle','on')
hold on
grid;
xlabel('Time (s)','FontSize',fs)
ylabel('lCM (molecules)','FontSize',fs)
plot(t,lCM,'-','Color',colorVec(iiBlack,:),'linewidth',lw)
if(savePic==1)
    saveFig(figname);
end

%% Solve Equations
fprintf(fLog, '%s: Solving Helical Rod model.\r\n', datestr(now,0));

R = sqrt(a0*Al/pi*m); % Basically, Al is a factor to take into account non-uniform distribution of protofibrils in fiber

[ls,le,Ft,T11,tanAlphaS,tanAlphaE,L,lSt] = ContinuumModel(...
    fLog,rm,a0,k0,tanAlpha0,K,Kb,tolContinuum,R,l,phiS,nuType,trackProgress...
    );


%% Output Parameters and Useful Values
fprintf(fLog, '%s: Outputting parameters and useful values.\r\n', datestr(now,0));
PiolaStress_OutputParams = filename+"_OutputParams.m";
fLogParams = fopen(PiolaStress_OutputParams,'w+');
if fLogParams == -1
    error('Cannot open log file.');
end
fprintf(fLogParams, '%% Time printed: %s.\r\n', datestr(now,0));
fprintf(fLogParams, '%% Source File Name:\r\n');
fprintf(fLogParams, 'StressODEModelFName = ''%s'' ;\r\n',filename);
fprintf(fLogParams, '%% Input Parameters:\r\n');
fprintf(fLogParams, 'rm = %g ;%% m\r\n', rm);
fprintf(fLogParams, 'k0 = %g ;%% m^(-1)\r\n', k0);
fprintf(fLogParams, 'mu = %g ;%% Pa s\r\n', mu);
fprintf(fLogParams, 'l = %g ;%% m\r\n', l);
fprintf(fLogParams, 'a0 = %g ;%% m^2\r\n', a0);
fprintf(fLogParams, 'Kb = %g ;%% N m^2\r\n', Kb);
fprintf(fLogParams, 'K = %g ;%% Pa\r\n', K);
if nuType==1 %fixed number density
    fprintf(fLogParams, 'nuType = %g ;%% (fixed number density)\r\n', nuType);
elseif nuType==2 %fixed solid volume fraction in free fiber config
    fprintf(fLogParams, 'nuType = %g ;%% (fixed solid volume fraction in free fiber config)\r\n', nuType);
elseif nuType==3 %fixed solid volume fraction in network config
    fprintf(fLogParams, 'nuType = %g ;%% (fixed solid volume fraction in network config)\r\n', nuType);
else
    %error('Invalid nuType %s.',nuType);
    fprintf(fLogParams, 'nuType = %g ;%% (unknown... script should have errored out prior to this anyway...)\r\n', nuType);
end
fprintf(fLogParams, '%% Output Parameters:\r\n');
fprintf(fLogParams, 'm = %g ; %% Average number of fiprotofibrils per fiber cros-section from Polymerization Model\r\n',m(end));
fprintf(fLogParams, 'n = %g ; %% Average number of fibrinogen monomers (fibrin) per protofibril from Polymerization Model (probably not correctly calculated)\r\n',n(end));
fprintf(fLogParams, 'lCM = %g ; %% Average length of fibers (in molecules) from Polymerization Model (probably not correctly calculated)\r\n',lCM(end));
fprintf(fLogParams, 'RE = %g ;%% m\r\n', R(end));
fprintf(fLogParams, 'L = %g ;%% m (=l/(lambdaE*lambdaS))\r\n', L);
fprintf(fLogParams, 'tanAlphaSn = %g ; \r\n', tanAlphaS(end));
fprintf(fLogParams, 'Fi = %g ;%% N\r\n', Ft(end));
fprintf(fLogParams, 'tanAlphaE = %g ;%% \r\n', tanAlphaE(end));
fprintf(fLogParams, 'lambdaS = %g ;%% \r\n', ls(end));
fprintf(fLogParams, 'lambdaE = %g ;%% \r\n', le(end));
fprintf(fLogParams, 'lambdaStar = %g ;%% \r\n', lSt(end));
fprintf(fLogParams, 'lend = %g ;%% m (=lambdaS*L)\r\n', ls(end)*L);
fprintf(fLogParams, 'TRzzE = %g ;%% Pa\r\n', T11(end));

%% Plotting Results
fprintf(fLog, '%s: Plotting results.\r\n', datestr(now,0));

fignum=fignum+1;
figname = filename+"__Figure_"+fignum; 
figure('Name',figname,'NumberTitle','on')
hold on
grid;
colorVector;
xlabel('\lambda^s','FontSize',fs)
ylabel('\lambda_2','FontSize',fs)
plot(ls(2:end),ls(2:end).*lSt(2:end),'-','Color',colorVec(iiBlack,:),'linewidth',lw)
if(savePic==1)
    saveFig(figname);
end

fignum=fignum+1;
figname = filename+"__Figure_"+fignum; 
figure('Name',figname,'NumberTitle','on')
hold on
grid;
colorVector;
xlabel('\lambda^s','FontSize',fs)
ylabel('T_{Rzz}','FontSize',fs)
plot(ls(2:end),T11(2:end),'-','Color',colorVec(iiBlack,:),'linewidth',lw)
if(savePic==1)
    saveFig(figname);
end

fignum=fignum+1;
figname = filename+"__Figure_"+fignum; 
figure('Name',figname,'NumberTitle','on')
hold on
grid;
colorVector;
xlabel('Time','FontSize',fs)
ylabel('T_{Rzz}','FontSize',fs)
plot(t(2:end),T11(2:end),'-','Color',colorVec(iiBlack,:),'linewidth',lw)
if(savePic==1)
    saveFig(figname);
end

if(haveExperimentalData==1) 
    fignum=fignum+1;
    figname = filename+"__Figure_"+fignum; 
    figure('Name',figname,'NumberTitle','on')
    hold on
    grid;
    colorVector;
    xlabel('Time','FontSize',fs)
    ylabel('T_{Rzz}','FontSize',fs)
    plot(t(2:end),T11(2:end),'-','Color',colorVec(iiBlack,:),'linewidth',lw)
    [AverageClotsWOSubtractionTime,AverageClotsWOSubtractionTension,...
        AverageControlTime,AverageClotsMinusControlTime,...
        AverageControlTension,AverageClotsMinusControlTension] ...
        =Import2021_07_14ExperimentalDataSets2();
    plot(AverageClotsWOSubtractionTime,AverageClotsWOSubtractionTension,'o','Color',colorVec(1,:),'linewidth',lw)
    plot(AverageControlTime,AverageControlTension,'o','Color',colorVec(2,:),'linewidth',lw)
    plot(AverageClotsMinusControlTime,AverageClotsMinusControlTension,'o','Color',colorVec(3,:),'linewidth',lw)
    legend({'Simulation','Average (Clots) Without Subtraction','Average Control (PPP, no clot)','Average (Clots) after Subtraction of Average Control'},'interpreter','latex','location','southeast')
    if(savePic==1)
        saveFig(figname);
    end
end

fignum=fignum+1;
figname = filename+"__Figure_"+fignum; 
figure('Name',figname,'NumberTitle','on')
hold on
grid;
colorVector;
xlabel('Time (s)','FontSize',fs)
ylabel('T_{Rzz} (Pa)','FontSize',fs)
plot(t(2:end),T11(2:end),'-','Color',colorVec(iiBlack,:),'linewidth',lw)
if(savePic==1)
    saveFig(figname);
end

fignum=fignum+1;
figname = filename+"__Figure_"+fignum; 
figure('Name',figname,'NumberTitle','on')
hold on
grid;
colorVector;
xlabel('Time','FontSize',fs)
ylabel('R (nm)','FontSize',fs)
plot(t(2:end),10^(9)*R(2:end),'-','Color',colorVec(iiBlack,:),'linewidth',lw)
if(savePic==1)
    saveFig(figname);
end

fignum=fignum+1;
figname = filename+"__Figure_"+fignum; 
figure('Name',figname,'NumberTitle','on')
hold on
grid;
xlabel('Time (s)','FontSize',fs)
ylabel('F (N)','FontSize',fs)
plot(t(2:end),Ft(2:end),'-','Color',colorVec(iiBlack,:),'linewidth',lw)
if(savePic==1)
    saveFig(figname);
end




%% Check PDE for Equilibrium Segment Length
if (doPDESegmentLengthEquilibrium==1)
    fprintf(fLog, '%s: Checking time to equilibrium Segment Length.\r\n', datestr(now,0));
    callingmfilename = filename;

    [tFiniteDifference,zl,fEndEl,s,ft0,ftPeak,ftOther,ftEnd,...
        res,restartwkspacename,MainWkspaceName,numRestarts...
    ] = RelaxationModel(...
        fLog,...
        PiolaStress_OutputParams,...
        callingmfilename,...
        tMin,tFinal,tRestartSize,...
        Dt,lenDiffEnd,...
        scriptfignum,savePic,trackProgress...
        restart,maxRestarts,restartWkspaceName,MainWkspaceName,maxArraySize...
        );
    
    matObjMain=matfile(char(MainWkspaceName),'Writable',true);
    tFiniteDifferenceEnd=eval(['matObjMain.t' num2str(numRestarts) '(1,end);']);
    zlEnd=eval(['matObjMain.zl' num2str(numRestarts) '(1,end);']);
    fprintf(fLogParams, '\r\n\r\ntFiniteDifferenceEnd = %g ;%% s\r\n',tFiniteDifferenceEnd);
    fprintf(fLogParams, 'zlEnd = %g ;%% s\r\n', zlEnd);
end




%% Save Workspace
wkspacename = filename+"_Workspace.mat";
save(wkspacename); % saves the workspace. can use whos('-file',wkspacename) to view contents, or load(wkspacename) to reload the workspace

%% Finish Up
tTotalFinish = toc(tTotalStart);
display([ ' Total script runtime: ' num2str(tTotalFinish)]);
fprintf(fLog, '%s: Total runtime: %f s.\r\n', datestr(now,0),tTotalFinish);
fprintf(fLog, '%s: Done.\r\n', datestr(now,0));
fclose(fLogParams);
fclose(fLog);
end