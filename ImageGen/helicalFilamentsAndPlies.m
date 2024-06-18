dbstop if error

tTotalStart=tic;
savePic = 1;
%doChecks = 0;

thisfilename=mfilename;

fLogThis = fopen(thisfilename+"_logFile.txt",'w+');
if fLogThis == -1
    error('Cannot open log file.');
end
fprintf(fLogThis, '%s: Loading input parameters for base ODE system.\r\n', datestr(now,0));

%color1=[0.4940 0.1840 0.5560];%purple
%{
color1=[0.4660 0.6740 0.1880];%green
color2=[0.8500 0.3250 0.0980];%orange
colorBase=[0.3010 0.7450 0.9330];%cyan
%}
% Chosen with triadic selector from https://www.sessions.edu/color-calculator/
% https://www.sessions.edu/color-calculator-results/?colors=9eddff,beee4d,ee4dbe
color1=[190 238 77]/255;
color2=[238 77 190]/255;
colorBase=[158 221 255]/255;%cyan

colors1=[[0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880]];
colors2=[[0 0.4470 0.7410];[0 0.4470 0.7410];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880]
    ];
colors3=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];
    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880]
    ];
colors4=[color1;color2;colorBase;colorBase;
    colorBase;colorBase;colorBase;colorBase;
    colorBase;colorBase;colorBase;colorBase;
    colorBase;colorBase;colorBase;colorBase;
    colorBase;colorBase;colorBase;colorBase;
    colorBase;colorBase;colorBase;colorBase;
    colorBase;colorBase;colorBase;colorBase;
    colorBase;colorBase;colorBase;colorBase;
    ];

rf=5;%nm
pf=400;%nm
rp=50;%nm
pp=1930;%nm

sfSize=2001;
sf=linspace(0,2000,sfSize);
theta=linspace(0,2*pi,32);
[sf,theta]=meshgrid(sf,theta);

nf=2;
np=16;

iff=1;
ip=1;

colors=colors4;
fignum=1;
%fignum=fignum+1;
figname = thisfilename+"__Figure_"+fignum;
figure('Name',figname,'NumberTitle','on','Units','normalized','Position',[0 0 1 1])
hold on
%grid;
%xlabel('Time (s)','FontSize',fs)
%ylabel('Stress (Pa)','FontSize',fs)
%plot(vals{1}{2}{2},vals{1}{2}{12},'-','Color',colorVec(iiBlack,:),'linewidth',lw)
%ii=1;plot(vals{2}{2}{2},vals{2}{2}{12},'-','Color',colorVec(ii,:),'linewidth',lw)
%ii=ii+1;plot(vals{3}{2}{2},vals{3}{2}{12},'-','Color',colorVec(ii,:),'linewidth',lw)
%legend({leg1,leg2,leg3,leg4,leg5},'interpreter','latex','location','southeast')
light;%('Position',[],'Style','infinite');
for ip=1:np
    for iff=1:nf
        %{
        surfaceCurve=zeros(3,len(sf),len(theta));
        for ii=1:len(theta)
            thetai=theta(ii);
            surfaceCurve(:,:,ii)=helicalFilamentsAndPliesFunction1(rf,pf,iff,nf,sf,rp,pp,ip,np,sf,rf,thetai);
        end
        %}
        surfaceCurve=helicalFilamentsAndPliesFunction1(rf,pf,iff,nf,sf,rp,pp,ip,np,sf,rf,theta);
        %{
        x=surfaceCurve(1,:,:);
        y=surfaceCurve(2,:);
        z=surfaceCurve(3,:);
        %}
        x=surfaceCurve(:,1:sfSize);
        y=surfaceCurve(:,sfSize+1:2*sfSize);
        z=surfaceCurve(:,2*sfSize+1:3*sfSize);
        surf(z,x,y,'FaceLighting','gouraud','EdgeColor','none','FaceColor',colors(2*ip-mod(iff,2),:));%'AmbientStrength',0.5?
    end
end
daspect([1 1 1]);
view(3);
set(gca,'visible','off');
if(savePic==1)
    saveFig(figname);
end


%% Save Workspace
filenameTogether = thisfilename;
wkspacename = filenameTogether+"_Workspace.mat";
save(wkspacename); % saves the workspace. can use whos('-file',wkspacename) to view contents, or load(wkspacename) to reload the workspace

%% Finish up
fprintf(fLogThis, '%s: Done.\r\n', datestr(now,0));
fclose(fLogThis);