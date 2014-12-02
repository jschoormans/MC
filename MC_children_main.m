%MONTE CARLO WITH MATING- MAIN FILE
clear all; close all; clc;

%GENERAL PARAMETERS
L=0.3;                              %arena size
Lb=0.01;                            %bot size
NbI=100;                            %initially created number of bots
Nb=10^3;                            %maximum number of bots;
Rm=0.038;                           %mating radius
t=100;                              %maximum simulation time
timeresolution=20;                  %number of timesteps per seconds=
lifetime=2.5;                       %lifetime of each bot
timespread=1;                       %initial time spread                   
survivalfraction=0.9;               %fraction of children that survives

%MATING PARAMETERS
infertility=1;                      %infertility switch
infertilityfrac=0.1;                %fraction of lifetime of infertility after mating;
maturity=1;                         %maturity switch
maturityfrac=0.25;                  %fraction of lifetime until maturity;
densitycheck=1;                     %switch for the density check mating requirements
densityradius=0.05;                 %radius of density check
densitynumber=25;                   %max number of other bots present in this radius to mate
speedfromparents=1;                 %if speed of bots is dependent on speed of  parents;


disp('Monte Carlo simulation for EvoComp - simulating bot mating')
disp('press any key to start')
pause;
%%


tic;
disp('Now running the simulation...')
[Activebots,BotSpeed,PosX,PosY,BotStartingTime]=MonteCarloBots(L,Lb,NbI,Nb,Rm,t,timeresolution,lifetime,timespread,survivalfraction,infertility,infertilityfrac,maturity,maturityfrac,densitycheck,densityradius,densitynumber,speedfromparents);
disp('Simulation finished')
toc

%% plot number of bots at any given time
figure(1);
plot(sum(Activebots))
xlabel('time')
ylabel('number of bots')

disp('Press any key to show bots')
% 
%% plot 
axis([0,0.5,0,0.5]);    

for ti=6:2:size(Activebots')
    close all

    figure(2);
    hold on
    plot(PosX(Activebots(:,ti),ti-5:ti)',PosY(Activebots(:,ti),ti-5:ti)','b.')  %plot part of bot traces
    scatter(PosX(:,ti)',PosY(:,ti)','r+')                                       %plot current bot position
    scatter(PosX(((BotStartingTime)==ti+1),ti+1)',PosY(((BotStartingTime)==ti+1),ti+1)','gp')
    %scatter(PosX(NbI+1:end,ti)',PosY(NbI+1:end,ti)','r+')
    text(0.28,0.28,strcat('t= ',num2str(ti/timeresolution),' '))
    
    hold off
   % plot([PosX(Parents(101,1),ti),PosX(Parents(101,2),ti)],[PosY(Parents(101,1),ti),PosY(Parents(101,2),ti)],'k')
   axis([0,0.3,0,0.3]);    
    drawnow
end


