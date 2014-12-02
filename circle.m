function [xoutput, youtput] = circle(x,y,r,BotAngSpeed,timesteps,BotOrient)


th = BotAngSpeed*timesteps+BotOrient;

xunit = r * cos(th)  ;

yunit = r * sin(th);

xoutput=xunit-xunit(1)+x;
youtput=yunit-yunit(1)+y;
