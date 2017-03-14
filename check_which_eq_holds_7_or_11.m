clear all; close all; clc

%calculate fmin and fplus for x direction, to check which equation applies
Ext = 658;
Rfitx = 514;
eta = 9.6E-10;
MLforcex = 4.9223;

Cpar = (1-9/16*(1+(Ext)/Rfitx)^(-1)+1/8*(1+(Ext)/Rfitx)^(-3)-45/256*(1+(Ext)/Rfitx)^(-4)-1/16*(1+(Ext)/Rfitx)^(-5))^(-1); %Daldrop eq(S10)
Crot = 1 + 5/16*(1+(Ext)/Rfitx)^(-3);
alphaX = 6*pi*eta*Rfitx*Cpar;% + 8*pi*eta*Rfitx*Crot/(1+Ext/Rfitx)^2; %Daldrop eq(11)
alphaPhi = 8*pi*eta*Rfitx^3*Crot; %Daldrop eq(11)

fPlus = (MLforcex/(Ext)*((Ext+Rfitx)*Rfitx/(2*alphaPhi) + 1/(2*alphaX) + 1/2*(((Ext+Rfitx)*Rfitx/alphaPhi + 1/alphaX)^2-4*Ext*Rfitx/(alphaX*alphaPhi))^(1/2)))/(2*pi);
fMin = (MLforcex/(Ext)*((Ext+Rfitx)*Rfitx/(2*alphaPhi) + 1/(2*alphaX) - 1/2*(((Ext+Rfitx)*Rfitx/alphaPhi + 1/alphaX)^2-4*Ext*Rfitx/(alphaX*alphaPhi))^(1/2)))/(2*pi);
