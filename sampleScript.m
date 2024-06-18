dbstop if error

filename = mfilename;
thisTic=tic;
fLog = fopen(filename+"_logFile.txt",'w+');
if fLog == -1
    error('Cannot open log file.');
end
fprintf(fLog, '%s: Starting %s.\r\n', datestr(now,0),filename);

%%
savefiles=1;
savePic=savefiles;
fignum=1;
lw=2;
fs=16;
linestyles={'-','-.','--',':'};

haveExperimentalData=0;
doPDESegmentLengthEquilibrium=1;
trackProgress=1;

%% Values for Piola Stress Calculations

% Constants in f
K=1.5*1000-4*139/3;%1.31467*10^3; %Pa %Bulk modulus of fibers, from Punter Soft Matter 2020

% Constants in eq90 (TRzz)
% L is given in the imaginary reference configuration.
% For the values in the network configuration [0.5, 0.75, 1, 2]*10^(-6),
% use the corresponding values of L in [0.8119, 1.2178, 1.6237, 2.4356, 3.2474]*10^(-6).
%L=0.8119*10^(-6); %m Length of fiber between crosslinks, in the imaginary reference configuration
% l is in the network configuration
l=0.5*10^(-6); %m Length of fiber between crosslinks, in the network configuration
%nu = 5*10^(15); %m^(-3) %3*sqrt(3)/L^3 from the 8-chain model. 
nuType=3;
% nuType 1: nu = 3*sqrt(3)/L^3; % No: of fibers per vol.
% nuType 2: nu = phiS./(pi*R^2*lambda_s*L); % (solid volume fraction in free fiber config)
% nuType 3: nu = phiS./(pi*R^2*lambda_s*lambda_e*L); % (solid volume fraction in network config)
phiS=0.01;
tolContinuum = 10^(-6); % tolerance of whether to accept solutions

% Constants in R
rm=6.5*10^(-9); %m Radius of fibrinogen monomers
a0=pi*(rm)^2; %m^2
% A crude method for introducing non-uniform distribution
% of protofibrils per fiber, as a fraction of m.
% A more rigorous method is described in the paper.
Al = 1; % Uniform distribution
%Al = 5/(13*sqrt(2));
%Al = (50*10^(-9)/rm)^2/200; 

% Constants in F
Kb=1600*10^(-12)*(10^(-9))^2; %N m^2 Bending modulus of fibers
k0=1.23*10^(-3)*1/10^(-9); %m^(-1)Spontaneous bending curvature of the helix
%alpha=82; %degrees
%alpha=deg2rad(alpha);

% Constants in lambdaS
tanAlpha0=400/(2*pi*5);

% Constants in Equilibrium Segment Length time sequence
mu = 1.002e-6*10^3;% Pa s. Viscocity of water

%% Values for Fiber Relaxation PDE Finite Difference Method

% Dt is an override parameter for default time step. To use the default, 
% Dt should be set to an empty character vector. To override the default,
% set Dt to whatever the override value should be.
%Dt =((6.238-6.162)/12)*10^(-08);
Dt = '';
tMin = 0.05*10^(-3);
tFinal = 1.0*10^(-3); % ms
tRestartSize = ''
lenDiffEnd = 0.01;% fraction difference from expected relaxed length at which to stop
scriptfignum=1;
restart=0;
maxRestarts=30;
restartWkspaceName='';
MainWkspaceName='';
maxArraySize=ceil(10^(-3)/(1.5262*10^(-10)));

%% Values for Polymerization ODE Set

xlimits = "auto"; %limits in the form [lower upper], or "auto"
ylimits = "auto"; %limits in the form [lower upper], or "auto"

tstep = 1;
tlim = 1000; %s
f_A0 = 2.8229; %mg/mL
% In the 1992 paper, it says a concentration of 0.5-1.5 molecule/L;
% Wikipedia (https://en.wikipedia.org/wiki/Fibrinogen) says the molecular 
% weight of fibrinogen is about 340 kDa, and also it 
% (https://en.wikipedia.org/wiki/Dalton_(unit)) says 
% 1 Da =  1.66053906660(50)*10^-27 kg.
% So I can convert 1 mg/mL= g/L 
%  = g/L (kg/1000g)(Da/1.66*10^-27 kg)(kDa/1000Da)(molecule/340 kDa), or 
% 1 mg/mL = 1/(1000*1000*340*1.66053906660(50)*10^(-27)) molecule/L
fibrinogenConversionFactor = 1/(1000*1000*340*1.66053906660*10^(-27));
%f_A0 = f_A0*fibrinogenConversionFactor; % Converts it from mg/mL to molecule/L
%f0 = 0;
f_n0 = 0;
f_r0 = 0;
f_nTot0 = 0;
k_A  = 1; %1/s
% Rule: k_pg > k_pi
k_pi = 6.0*10^(-20); %L/molecule s (multiply by Avogadro's number to get moles?)
k_pg = 1.4*10^(-17); %L/molecule s (multiply by Avogadro's number to get moles?)
% Rule: k_fg > k_fi
k_fi = 1*10^(-20); %L/molecule s (multiply by Avogadro's number to get moles?)
k_fg = 2*10^(-16); %L/molecule s (multiply by Avogadro's number to get moles?)
k_fA = 1*10^(-19);
c_fn0 = 0;
c_fr0 = 0;
minLength = 20;
f0 = zeros(1,minLength);


%% Changes for variable, and computation loop:

param='$l$';
%paramVar='double(subs(lambda_ce*ls(end)*L,[lambdaStar,lambda_s,Nm],[lSt(end),lsR(end),rangeParam(end)]))';
paramVarName='l';
unit='$m$';
base=filename+'_';
p = [2,1.5,1,0.75,0.5].*10^(-6);
params=cell(size(p));
params{1}={p{1},[base paramVarName '200']};
params{2}={p{2},[base paramVarName '150']};
params{3}={p{3},[base paramVarName '100']};
params{4}={p{4},[base paramVarName '075']};
params{5}={p{5},[base paramVarName '050']};
outputs=cell(size(p));

fLogMain=fLog;

for ii=1:length(params)
    fprintf(fLogMain, '%s: Running %s.\r\n', datestr(now,0),params{ii}{2});
    disp(['Running ' params{ii}{2}]);
    l=params{ii}{1};
    loopfilename=params{ii}{2};
    wkspacename = PiolaStressFromPolymerization(...
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
        loopfilename...
        );
    load(wkspacename,'-mat')
    %l=eval(paramVar);
    if(doPDESegmentLengthEquilibrium==1)
        outputs{ii}={...
            t,m,n,lCM,...
            T11,...
            tFiniteDifference,zl,fEndEl,s,ft0,ftPeak,ftOther,ftEnd...
        }
    else
        outputs{ii}={...
            t,m,n,lCM,...
            T11...
        }
    end
end



%% Save Workspace
wkspacename = filename+"_Workspace.mat";
save(wkspacename); % saves the workspace. can use whos('-file',wkspacename) to view contents, or load(wkspacename) to reload the workspace



%% Plot Results
fprintf(fLogMain, '%s: Plotting results together.\r\n', datestr(now,0));

leg=cell(size(p));
colorVector;

fignum=fignum+1;
figname = filename+"__Figure_"+fignum; 
figure('Name',figname,'NumberTitle','on')
hold on
grid;
xlabel('Time (s)','FontSize',fs)
ylabel('$T_{Rzz}~Pa$','FontSize',fs,'interpreter','latex')
%title('$T_{Rzz} vs. t$','interpreter','latex','FontSize',fs)
for ii=1:length(p)
    plot(outputs{ii}{1},outputs{ii}{5},'-','Color',colorVec(ii,:),'linewidth',lw,'linestyle',linestyles{ii})
    leg{ii}=char(sprintf('%s $=%0.2g$ %s',param,params{ii}{1},unit));
end
legend(leg,'interpreter','latex','Position',[0.6301 0.5532 0.2515 0.1857]);
if(savePic==1)
    saveFig(figname);
end


if(doPDESegmentLengthEquilibrium==1)
    fignum=fignum+1;
    figname = mainfilename+"__Figure_"+fignum; 
    figure('Name',figname,'NumberTitle','on')
    hold on
    grid;
    xlabel('Time (s)','FontSize',fs)
    ylabel('$z(s=l,t)$','interpreter','latex','FontSize',fs)
    for ii=1:length(p)
        plot(outputs{ii}{6},outputs{ii}{7},'-','Color',colorVec(ii,:),'linewidth',lw,'linestyle',linestyles{ii})
        leg{ii}=char(sprintf('%s $=%0.2g$ %s',param,params{ii}{1},unit));
    end
    legend(leg,'interpreter','latex','location','southeast');
    if(savePic==1)
        saveFig(figname);
    end

    fignum=fignum+1;
    figname = mainfilename+"__Figure_"+fignum; 
    figure('Name',figname,'NumberTitle','on')
    hold on
    grid;
    xlabel('Time (s)','FontSize',fs)
    ylabel('$\frac{z(s=l,t)}{z(s=l,t=0)}$','interpreter','latex','FontSize',fs)
    for ii=1:length(p)
        plot(outputs{ii}{6},outputs{ii}{7}./outputs{ii}{7}(1),'-','Color',colorVec(ii,:),'linewidth',lw,'linestyle',linestyles{ii})
        leg{ii}=char(sprintf('%s $=%0.2g$ %s',param,params{ii}{1},unit));
    end
    legend(leg,'interpreter','latex','location','southeast');
    if(savePic==1)
        saveFig(figname);
    end
end


%% Finish up
thisToc=toc(thisTic);
message = sprintf('Total script time: %s s.',num2str(thisToc));
fprintf('%s: Done with %s. %s\r\n', datestr(now,0),filename,message);
fprintf(fLogMain, '%s: Done with %s. %s\r\n', datestr(now,0),filename,message);
fclose(fLogMain);