function [t,m,n,l,varargout] = PolymerizationModel(params)

fname = params{1};
savePic = params{2};

fLog = fopen(fname+"logFile.txt",'w+');
if fLog == -1
    error('Cannot open log file.');
end
%% Input parameters
fprintf(fLog, '%s: Loading input parameters.\r\n', datestr(now,0));

tstep = params{3};
tlim = params{4};
f_A0 = params{5};
f0 = params{6};
f_n0 = params{7};
f_r0 = params{8};
f_nTot0 = params{9};
k_A  = params{10};
k_pi = params{11};
k_pg = params{12};
k_fi = params{13};
k_fg = params{14};
c_fn0 = params{17};
c_fr0 = params{18};
minLength = params{19};
k_fA = params{20};
if minLength<2
    minLength=2;
    warning('minLength must be >= 2. Setting minLength=2.');
end

%% ODE parameters
fprintf(fLog, '%s: Setting up ODE parameters.\r\n', datestr(now,0));

tspan = [0.0 : tstep : tlim]';
yDe0 = zeros(1,minLength+6);
yDe0(1) = f_A0;
%yDe0(2) = f0;
for ii=2:minLength+1
    yDe0(ii) = f0(ii-1);
end
yDe0(minLength+2) = f_n0;
yDe0(minLength+3) = f_r0;
yDe0(minLength+4) = f_nTot0;
yDe0(minLength+5) = c_fn0;
yDe0(minLength+6) = c_fr0; 

dotParams{1} = fLog;
dotParams{2} = k_A;
dotParams{3} = k_pi;
dotParams{4} = k_pg; 
dotParams{5} = k_fi;
dotParams{6} = k_fg;
dotParams{7} = minLength;
dotParams{8} = k_fA;

%% Solving ODE 
fprintf(fLog, '%s: Solving ODE system.\r\n', datestr(now,0));

tic;
clear textprogressbar
ops = odeset('OutputFcn',@odePolymerizationModeltpbar);
[t,yDe] = ode45(@PolymerizationModelOdeSet,tspan,yDe0,ops,dotParams);

f_A    = yDe(:,1);
f      = yDe(:,2:minLength+1);
f_n    = yDe(:,minLength+2);
f_r    = yDe(:,minLength+3);
f_nTot = yDe(:,minLength+4);
c_fn   = yDe(:,minLength+5);
c_fr   = yDe(:,minLength+6);
m = f_nTot./f_r; % Average number protofibril/fiber cross-section
n = c_fn./f_n;   % Average number fibrin/protofibril
l = c_fr./f_nTot;% Average length of fibers

m(1) = 0;
n(1) = 0;
l(1) = 0;

nout = max(nargout,4)-4;
if nout>0
    if nout==1
        varargout={f_A,f,f_n,f_r,f_nTot,c_fn,c_fr};
    else
        varargout{1}=f_A;
        varargout{2}=f;
        if nout>2
            varargout{3}=f_n;
        end
        if nout>3
            varargout{4}=f_r;
        end
        if nout>4
            varargout{5}=f_nTot;
        end
        if nout>5
            varargout{6}=c_fn;
        end
        if nout>6
            varargout{7}=c_fr;
        end
    end
end
   

%% Plotting results
fprintf(fLog, '%s: Plotting results.\r\n', datestr(now,0));

lw=2;
fs = 16;
figure('Name',fname,'NumberTitle','on')
hold on
colorVector;
ii=1;
xlabel('time (s)','FontSize',fs)
ylabel('protofibrils per fiber','FontSize',fs)

plot(t,m,'-','Color',colorVec(iiBlack,:),'linewidth',lw)

%if params{15} ~= "auto"
if ~isstring(params{15})
    if isvector(params{15})
        xlim(params{15});
    end
end
%if params{16} ~= "auto"
if ~isstring(params{16})
    if isvector(params{16})
        ylim(params{16});
    end
end

if(savePic==1)
    fprintf(fLog, '%s: Saving files.\r\n', datestr(now,0));
    saveFig(fname);
end

%% Finish
fprintf(fLog, '%s: Done.\r\n', datestr(now,0));
fclose(fLog);
end