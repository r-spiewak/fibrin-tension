function [t,zl,fEndEl,s,ft0,ftPeak,ftOther,ftEnd,res,...
    restartwkspacename,MainWkspaceName,numRestarts] = RelaxationModel(...
        fLog,...
        PiolaStress_OutputParams,...
        callingmfilename,...
        tMin,tFinal,tRestartSize,...
        Dt,lenDiffEnd,...
        scriptfignum,savePic,trackProgress,...
        restart,maxRestarts,restartWkspaceName,MainWkspaceName,maxArraySize...
        )
dbstop if error
makeRestart=0;
restartCreated=0;
thisFuncName=mfilename;
doPlotting=0;
plotted=0;
clear textprogressbar

run(char(PiolaStress_OutputParams));

if ~restart
    Ds = 1*10^(-9);
    s=linspace(0,l,ceil(l/Ds));
    % Increments of time must be smaller than this:
    % (s(2)-s(1))/(2*2*pi*mu/(log(L/R)+log(2)-3/2))
    d_w = double(2*pi*mu/(log(l/RE)+log(2)-3/2));
    Dt1 = Dt;
    hypS = sqrt(1+tanAlphaSn.^2);
    sinAlphaS = tanAlphaSn./hypS;
    cosAlphaS = 1./hypS;
    Ka = 2*pi/a0*Kb*tanAlphaE*cosAlphaS^3*(k0*(RE-rm)+(1+sinAlphaS^2)*log(RE/rm));
    Dt2 = d_w/(2*Ka)*Ds^2;
    %if Dt2<Dt1 
    % In theory, this is the correct condition for convergence. However, 
    % it does not seem to depend on
    % anything that differs for different L, and I remember that some of the
    % bigger L had issues with this Dt (but I cannot find proof of that, at 
    % this moment), so this may not work for everything. But I will try it
    % anyway and see. I can always modify it again later...
    %if Dt2>Dt1
    % Somehow there is dependence, even if I don't see it, so I will just
    % change the condition to use Dt1 only if the input isn't blank (or
    % something like that).
    if isempty(Dt1)
        Dt=Dt2;
    end
    try
        t=linspace(0,tFinal,ceil(tFinal/Dt));
        tRestartSize=tFinal;
        tFinalNew=tFinal;
    catch ME1
        idSegLast = regexp(ME1.identifier,'(?<=:)\w+$','match');
        %keyboard
        if strcmp(idSegLast,'pmaxsize')||strcmp(idSegLast,'SizeLimitExceeded')
            message = sprintf('%s\nWill instead split time interval and \ncreate Restart file if necessary.',ME1.message);
            warning(message)
            fprintf(fLog, '%s: Warning: %s.\r\n', datestr(now,0),message);
            tRestartSize=tFinal/2;
            while (ceil(tRestartSize/Dt)>maxArraySize)
                tRestartSize=tRestartSize/2;
            end
            tFinalNew=tRestartSize;
            t=linspace(0,tRestartSize,ceil(tRestartSize/Dt));
        else
            keyboard
            rethrow(ME1);
        end
    end
    MainWkspaceName = callingmfilename+"_MainWorkspace.mat";
    save(char(MainWkspaceName),'-v7.3');
    matObjMain = matfile(char(MainWkspaceName),'Writable',true);
else
    fLogSaved = fLog;
    callingmfilenameSaved = callingmfilename;
    %load some things from Main workspace
    matObjMain = matfile(char(MainWkspaceName),'Writable',true);
    fLog=matObjMain.fLog;
    callingmfilename=matObjMain.callingmfilename;
    tFinal=matObjMain.tFinal;
    Dt=matObjMain.Dt;
    s=matObjMain.s;
    Ds=matObjMain.Ds;
    zOld=matObjMain.zOld;
    d_w=matObjMain.d_w;
    tRestartSize=matObjMain.tRestartSize;
    if ~isequal(fLogSaved,fLog)
        fLog=fLogSaved;
    end
    if ~strcmp(callingmfilenameSaved,callingmfilename)
        callingmfilename=callingmfilenameSaved;
    end
    tRestart=eval(['matObjMain.t' num2str(restart-1) '(1,end);']);
    if tRestart== tFinal
        tFinal=tFinal+tRestartSize;
        tFinalNew=tFinal;
    else
        tFinalNew=tRestart+tRestartSize;
    end
    try
        t=linspace(tRestart,tFinalNew,ceil((tFinalNew-tRestart)/Dt));
    catch ME1
        idSegLast = regexp(ME1.identifier,'(?<=:)\w+$','match');
        keyboard
        if strcmp(idSegLast,'pmaxsize')||strcmp(idSegLast,'SizeLimitExceeded') 
            message = sprintf('%s\nWill instead split time interval and \ncreate Restart file if necessary.',ME1.message);
            warning(message)
            fprintf(fLog, '%s: Warning: %s.\r\n', datestr(now,0),message);
            %makeRestart=1;
            keyboard
            tFinalNew=tFinalNew/2;
            t=linspace(tRestart,tFinalNew,ceil((tFinalNew-tRestart)/Dt));
        else
            keyboard
            rethrow(ME1);
        end
    end
end
%% Here I implement the Finite Difference method
% s represents the nodes in my rod; elements exist between these nodes
%zinit = [0:Ds:l];
Omega1 = tanAlphaE;
% This is the simple initial condition, which isn't actually consistent
% with the boundary conditions. Presumably we should change it eventually.
zinit = s;
for j=1:length(s)
    if (s(j)>0.99*l)
        zinit(j) = (tanAlphaSn/Omega1)*(s(j)-0.99*l)+0.99*l;
    end
end
Ns = length(s); %number of nodes
ds = Ds;
% include the boundary conditions probably in the loop (except for t=0)
%z=zeros(length(t),length(zinit));
%z(1,:) = zinit;
zl = 0*t;
if ~restart
    zOld = zinit*1;
    zl(1) = zinit(end);
else
    zl(1) = zOld(end);
end
fEndEl = 0*t;
zlrelaxed = l*tanAlphaSn/Omega1;
itPeak = find(t>10^(-8),1,'first');
itOther = find(t>10^(-7),1,'first');
if(trackProgress==1)
    tstart=tic;
    tfin=length(t);
    fprintf('Restart: %d\n',restart);
    textprogressbar('Finite Difference Method progress: ');
end
for j=2:length(t) 
    % loop stuff here
    tanAlpha = Omega1 * (zOld(2:end)-zOld(1:end-1)) / ds; % tanA
    hyp = sqrt(1+tanAlpha.^2);
    sinAlpha = tanAlpha./hyp;
    cosAlpha = 1./hyp;
    F = (2*pi/a0)*Kb*sinAlpha.*(-log(RE./rm).*(cosAlpha).^2+k0.*(RE-rm));
    if j==2
        ft0 = F;
    elseif j-1==itPeak
        ftPeak = F;
    elseif j-1==itOther
        ftOther = F;
    end
    %F = ((2*pi/a0)*Kb)*((z(j-1,2:end)-z(j-1,1:end-1)) / ds -1); %Linear
    %force-stretch relation for testing purposes...
    % The above are the forces in the elements, or between two nodes.
    dt = t(j)-t(j-1);
    dz = dt*(F(2:end)-F(1:end-1))/(d_w*ds);
    dz1 = 0;   %% Node 1 is fixed
    dzend = (0-F(end))*dt/(d_w*ds);
    dzcur = [dz1 dz dzend]; 
    zOld = zOld + dzcur;
    fEndEl(j-1) = F(end);
    zl(j) = zOld(end);
    % Something here to be able to track progress of the scheme.
    if (trackProgress == 1)
        progress=100*j/tfin;
        textprogressbar(progress);
    end
    if abs(zlrelaxed-zl(j))/zlrelaxed<lenDiffEnd && t(j)>tMin
        break
    end
end
if (trackProgress==1)
    tfin=[];
    textprogressbar('');
    display([ '   Finite Difference Method time: ' num2str(toc(tstart))]);
    tstart=[];
end
tanAlpha = Omega1 * (zOld(2:end)-zOld(1:end-1)) / ds; % tanA
hyp = sqrt(1+tanAlpha.^2);
sinAlpha = tanAlpha./hyp;
cosAlpha = 1./hyp;
F = (2*pi/a0)*Kb*sinAlpha.*(-log(RE./rm).*(cosAlpha).^2+k0.*(RE-rm));
ftEnd = F;
fEndEl(j) = F(end);
%if restart
    if exist('ftPeak','var')==0
        ftPeak='';
    end
    if exist('ftOther','var')==0
        ftOther='';
    end
%end
if t(j)<t(end)
    t=t(1:j);
    zl=zl(1:j);
    fEndEl=fEndEl(1:j);
    numRestarts=restart;
    matObjMain.numRestarts=numRestarts;
    makeRestart=0;
else
    message = sprintf('PDE output length did not relax to less than tolerance %g in time %gs. \n zlRelaxed=%g, zl(j)=%g, abs(zlrelaxed-zl(j))/zlrelaxed=%g.',lenDiffEnd,tFinalNew,zlrelaxed,zl(j),abs(zlrelaxed-zl(j))/zlrelaxed);
    warning(message)
    fprintf(fLog, '%s: Warning: %s\r\n', datestr(now,0),message);
    makeRestart=1;
end
eval(['matObjMain.t' num2str(restart) '=t*1;']);
eval(['matObjMain.zl' num2str(restart) '=zl*1;']);
eval(['matObjMain.fEndEl' num2str(restart) '=fEndEl*1;']);
matObjMain.zlrelaxed=zlrelaxed;
matObjMain.fLog=fLog;
matObjMain.callingmfilename=callingmfilename;
matObjMain.tFinal=tFinal;
matObjMain.tRestartSize=tRestartSize;
matObjMain.maxArraySize=maxArraySize;
matObjMain.zOld=zOld;
if (makeRestart)
    restart=restart+1;
    restartName=callingmfilename+"_restart"+num2str(restart);
    message = sprintf('Creating Restart files. Find file %s',restartName);
    warning(message)
    fprintf(fLog, '%s: %s.\r\n', datestr(now,0),message);
    matObjMain.restart=restart;
    restartWkspaceName = restartName+"_Workspace.mat";
    MainFullFileName=fullfile(pwd,char(MainWkspaceName));
    restartFullFileName=fullfile(pwd,char(restartWkspaceName));
    message = sprintf('Backing up:\n %s \nto\n %s.\n',MainWkspaceName,restartWkspaceName);
    fprintf(fLog, '%s: %s.\r\n', datestr(now,0),message);
    fprintf('%s\n',message);
    copyfile(char(MainFullFileName),char(restartFullFileName),'f');
    restartCreated=1;
end
if (restartCreated)
    %if not run, inform (and log) how restart can be run later.
    gotoend=0;
    message=sprintf('Restart can be run with command following key:\n [t,zl,fEndEl,s,ft0,ftPeak,ftOther,ftEnd,res,restartwkspacename,MainWkspaceName,numRestarts]=%s(fLog,PiolaStress_OutputParams,callingmfilename,tMin,tFinal,tFinalNew,Dt,lenDiffEnd,scriptfignum,savePic,trackProgress,restart,maxRestarts,restartWkspaceName,MainWkspaceName,maxArraySize)\n [t,zl,fEndEl,s,ft0,ftPeak,ftOther,ftEnd,res,restartwkspacename,MainWkspaceName,numRestarts]=%s(%d,%s,%s,%g,%g,%g,%g,%f,%d,%d,%d,%d,%d,"%s","%s",%g)',thisFuncName,thisFuncName,fLog,PiolaStress_OutputParams,callingmfilename,tMin,tFinal,tRestartSize,Dt,lenDiffEnd,scriptfignum,savePic,trackProgress,restart,maxRestarts,restartWkspaceName,MainWkspaceName,maxArraySize);
    warning(message)
    fprintf(fLog, '%s: %s.\r\n', datestr(now,0),message);
    doPlotting=0;
    if restart<=maxRestarts
        %ask if it should run restart, and run it?
        question='Should the Restart be run automatically?';
        button=questdlg_timer(10,question,'Run Restart?');
        if strcmp(button,'Yes')
            message=sprintf('Running Restart with command following key:\n [t,zl,fEndEl,s,ft0,ftPeak,ftOther,ftEnd,res,restartwkspacename,MainWkspaceName,numRestarts]=%s(fLog,PiolaStress_OutputParams,callingmfilename,tMin,tFinal,tFinalNew,Dt,lenDiffEnd,scriptfignum,savePic,trackProgress,restart,maxRestarts,restartWkspaceName,MainWkspaceName,maxArraySize)\n [t,zl,fEndEl,s,ft0,ftPeak,ftOther,ftEnd,res,restartwkspacename,MainWkspaceName,numRestarts]=%s(%d,%s,%s,%g,%g,%g,%g,%f,%d,%d,%d,%d,%d,"%s","%s",%g)',thisFuncName,thisFuncName,fLog,PiolaStress_OutputParams,callingmfilename,tMin,tFinal,tRestartSize,Dt,lenDiffEnd,scriptfignum,savePic,trackProgress,restart,maxRestarts,restartWkspaceName,MainWkspaceName,maxArraySize);
            warning(message)
            fprintf(fLog, '%s: %s.\r\n', datestr(now,0),message);
            [t,zl,fEndEl,s,ft0,ftPeak,ftOther,ftEnd,res,restartwkspacename,MainWkspaceName,numRestarts]=eval([thisFuncName '(fLog,PiolaStress_OutputParams,callingmfilename,tMin,tFinal,tFinalNew,Dt,lenDiffEnd,scriptfignum,savePic,trackProgress,restart,maxRestarts,restartWkspaceName,MainWkspaceName,maxArraySize);']);
            doPlotting=0;
        else
            gotoend=1;
        end
    else
        %Inform user (i.e., me) that maxRestarts has been reached.
        message=sprintf('Maximum number of restarts have already been run. Increase maxRestarts to run more Restarts.');
        warning(message)
        fprintf(fLog, '%s: Warning: %s.\r\n', datestr(now,0),message);
        doPlotting=0;
        gotoend=1;
    end
    if gotoend==1
        % gotoend should presumably do something... but I guess it's not
        % necessary. I'll let it make this output variable though.
        numRestarts=restart;
        matObjMain.numRestarts=numRestarts;
    end
end


%% Plot figures
if doPlotting==1
    if exist('scriptfignum','var')==0
        scriptfignum=1;
    else
        scriptfignum=scriptfignum+1;
    end
    lw=2;
    fs = 16;
    wkspacename = callingmfilename+"_pdeFiniteDifference_Workspace.mat";
    save(wkspacename); % saves the workspace. can use whos('-file',wkspacename) to view contents, or load(wkspacename) to reload the workspace
    figname = callingmfilename+"__pdeFiniteDifference_Figure_"+scriptfignum; 
    figure('Name',figname,'NumberTitle','on')
    hold on;
    for ii=0:restart
        tArray=eval(['matObjMain.t' num2str(ii)]);
        zlArray=eval(['matObjMain.zl' num2str(ii)]);
        plot(tArray,zlArray,'b-o')
    end
    xlabel('$t$ (s)','FontSize',fs,'interpreter','latex')
    ylabel('$z(s=l,t)$ (m)','FontSize',fs,'interpreter','latex')
    if(savePic==1)
        saveFig(figname);
    end
    scriptfignum=scriptfignum+1;
    figname = callingmfilename+"__pdeFiniteDifference_Figure_"+scriptfignum; 
    figure('Name',figname,'NumberTitle','on')
    hold on;
    for ii=0:restart
        tArray=eval(['matObjMain.t' num2str(ii)]);
        fEndElArray=eval(['matObjMain.fEndEl' num2str(ii)]);
        plot(tArray,fEndElArray,'b-o')
    end
    xlabel('$t$ (s)','FontSize',fs,'interpreter','latex')
    ylabel('$F(s=l,t)$ (N)','FontSize',fs,'interpreter','latex')
    if(savePic==1)
        saveFig(figname);
    end
end
res = scriptfignum;
if exist('restartwkspacename','var')==0
    restartwkspacename=restartWkspaceName;
end
end
