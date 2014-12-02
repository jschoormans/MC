%%MONTE CARLO MATING WITH OFFSPRING!!!
clear all; close all; clc;

tic
L=0.3;                          %arena size
A=L*L;                          %Area of Arena
Lb=0.01;                        %bot size
NbI=100;                        %initially created number of bots
Nb=10^4;                        %maximum number of bots;
Rm=0.038;                       %mating radius
t=100;                         %maximum simulation time
timeresolution=20;             %number of timesteps per seconds=
NrTimeSteps=timeresolution*t;
lifetime=2.5;                   %lifetime of each bot
timespread=1;                   %initial time spread                   
timesteps=[1:NrTimeSteps];      %timestep vector
timestepsbot=[1:lifetime*timeresolution];
survivalfraction=0.9;           %fraction of children that survives
timevector=linspace(0,t,NrTimeSteps)

infertility=1;                  %infertility switch
infertilityfrac=0.1;           %fraction of lifetime of infertility after mating;
maturity=1;                     %maturity switch
maturityfrac=0.25               %fraction of lifetime until maturity;
densitycheck=1;
densityradius=0.05;
densitynumber=25;
speedfromparents=1;             %if speed of bots is dependent on speed of  parents;



%initalize Bot Posistions
BotPositionX=zeros(Nb,NrTimeSteps);     
BotPositionY=zeros(Nb,NrTimeSteps);
BotOrient=zeros(Nb,NrTimeSteps);    %orientation of Bot in radians
BotRadius=zeros(Nb,1);
BotSpeed=zeros(Nb,1);

BotStartingTime=randi(timespread*timeresolution,NbI,1); %starting time of initial bots



%non-dynamic bot-specific parameters (for all bots, also children!)

load('speedparameters')
BotSpeed=exp(erfinv(2*(rand(Nb,1)-0.5))*(sqrt(2)*lsigma)+lmu);
BotRadius=abs(randn(Nb,1)*(L/4)); %specify radii of trajectories
BotOrient=rand(Nb,1)*(2*pi);
BotAngSpeed=(BotSpeed./(2*pi*BotRadius))*(2*pi)*(t/NrTimeSteps);

%dynamic bot-specific parameters
for i=1:NbI
BotPositionX(i,BotStartingTime(i))=rand(1,1)*L;
BotPositionY(i,BotStartingTime(i))=rand(1,1)*L;

end



%% CALCULATE  POSITIONS FOR ALL TIMESTEPS
PosX=nan(Nb,NrTimeSteps);
PosY=nan(Nb,NrTimeSteps);
for n=1:NbI
[PosX(n,BotStartingTime(n):BotStartingTime(n)+(length(timestepsbot)-1)),PosY(n,BotStartingTime(n):BotStartingTime(n)+(length(timestepsbot)-1))]=circle(BotPositionX(n,BotStartingTime(n)),BotPositionY(n,BotStartingTime(n)),BotRadius(n),BotAngSpeed(n),timestepsbot,BotOrient(n,1));
end

%%
%BOUNDARY CONDS.
while sum(sum(PosX(1:NbI,:)>L))>0
    A=PosX(1:NbI,:)>L;
    PosX(1:NbI,:)=PosX(1:NbI,:)-A.*L;
end
while sum(sum(PosX(1:NbI,:)<0))>0;
    A=PosX(1:NbI,:)<0;
    PosX(1:NbI,:)=PosX(1:NbI,:)+A.*L;
end

while sum(sum(PosY(1:NbI,:)>L))>0
    A=PosY(1:NbI,:)>L;
    PosY(1:NbI,:)=PosY(1:NbI,:)-A.*L;
end
while sum(sum(PosY(1:NbI,:)<0))>0;
    A=PosY(1:NbI,:)<0;
    PosY(1:NbI,:)=PosY(1:NbI,:)+A.*L;
end



%% NOW: FIND MATING INSTANCES
NMatings=zeros(NrTimeSteps,1);    %number of mating events per timestep
Matings=diag(ones(Nb,1));           %matrix for parent control
NewBotId=NbI+1;                     %ID for newly created bot
survival=rand(Nb,1)<survivalfraction; %vector of surviving children
Parents=nan(Nb,2);
Matingtime=zeros(Nb,1);


tt=1;
while NewBotId<Nb & tt<NrTimeSteps;    
    
Activebots(:,tt)=isnan(PosX(:,tt))==0;
ActivebotsID=find(Activebots(:,tt)==1);
    if sum(Activebots(:,tt))>1;         %no distances can be calculated from just one or zero bots
        distances=(pdist([PosX(Activebots(:,tt),tt),PosY(Activebots(:,tt),tt)]));
        D=squareform(distances);        
        [D1,D2]=size(D);
         for i=2:D1;
             for j=1:i-1;

             BotId1=ActivebotsID(i);
             BotId2=ActivebotsID(j);
             
                 if (maturity==1)&((tt-BotStartingTime(BotId1))>(infertilityfrac*lifetime*timeresolution))&((tt-BotStartingTime(BotId2))>(infertilityfrac*lifetime*timeresolution)) || maturity==0;
                 if (infertility==1)& ((tt-Matingtime(BotId1))>(infertilityfrac*lifetime*timeresolution))&((tt-Matingtime(BotId2))>(infertilityfrac*lifetime*timeresolution)) || infertility==0;
                        
                     if D(i,j)<Rm && Matings(BotId1,BotId2)==0;
                         if densitycheck==1&sum(D(i,:)<densityradius)<densitynumber&sum(D(j,:)<densityradius)<densitynumber||densitycheck==0;
                     if NewBotId<Nb;
                     Matings(BotId1,BotId2)=1;  %mating constraints
                     Matings(BotId2,BotId1)=1;

                     Matings(NewBotId,BotId1)=1;
                     Matings(NewBotId,BotId2)=1;
                     Matings(BotId1,NewBotId)=1;
                     Matings(BotId2,NewBotId)=1;
                     
                     Matingtime(BotId1)=tt;     %save time of mating
                     Matingtime(BotId2)=tt;     %save time of mating

                     NMatings(tt)=NMatings(tt)+1;

                     %%ADD NEW BOT:
                     if survival(NewBotId)==1;
                        BotStartingTime(NewBotId)=tt+1;
                        BotStartingTime(NewBotId);
                        
                        Xstart=(PosX(BotId1,tt)+PosX(BotId2,tt))/2;
                        Ystart=(PosY(BotId1,tt)+PosY(BotId2,tt))/2;
                        Parents(NewBotId,:)=[BotId1,BotId2];
                        
                        if speedfromparents==1;     %Bot speed depends on speed of parents (to mimic evolution)
                        BotSpeed(NewBotId)=mean([BotSpeed(BotId1),BotSpeed(BotId2)])+randn*(max([BotSpeed(BotId1),BotSpeed(BotId2)])-mean([BotSpeed(BotId1),BotSpeed(BotId2)]));           %bot speed is average of parent's speed (NOT SURE IF GOOD TO USE?)
                        end
                        
                        [PosX(NewBotId,BotStartingTime(NewBotId):BotStartingTime(NewBotId)+(length(timestepsbot)-1)),PosY(NewBotId,BotStartingTime(NewBotId):BotStartingTime(NewBotId)+(length(timestepsbot)-1))]=circle(Xstart,Ystart,BotRadius(NewBotId),BotAngSpeed(NewBotId),timestepsbot,BotOrient(NewBotId,1));
                        %ADD BOUNDARY CONDS
                    
                        while (sum(PosX(NewBotId,:)>L))>0
                            A=PosX(NewBotId,:)>L;
                            PosX(NewBotId,:)=PosX(NewBotId,:)-A.*L;
                        end
                        while (sum(PosX(NewBotId,:)<0))>0;
                            A=PosX(NewBotId,:)<0;
                            PosX(NewBotId,:)=PosX(NewBotId,:)+A.*L;
                        end

                        while (sum(PosY(NewBotId,:)>L))>0
                            A=PosY(NewBotId,:)>L;
                            PosY(NewBotId,:)=PosY(NewBotId,:)-A.*L;
                        end
                        while (sum(PosY(NewBotId,:)<0))>0;
                            A=PosY(NewBotId,:)<0;
                            PosY(NewBotId,:)=PosY(NewBotId,:)+A.*L;
                        end
                     end
                     end
                    NewBotId=NewBotId+1;
                     end
                 end
             end
             end
             end
        end
clear distances, clear D;
    end    
    tt=tt+1;
end
toc
%% plot number of bots at any given time

figure(1);
plot(sum(Activebots))
xlabel('time')
ylabel('number of bots')

figure(2)
hold on
movavg(BotSpeed,1000,2000,1)
hold off
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
