%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Top level section 2: FX2_LOAD_RAW  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%% This is the second section in a series of routines to analyze
    %%% force-extension data taken with the new multi bead code for MT
    %%% (J. P. Cnossen, D. Dulin and N. H. Dekker Rev Sci Instrum 85, (2014).)
    %%%
    %%% This section is used to load the data for the force-extension
    %%% measurements, namely the large matrix of position traces for different
    %%% magnet heights.
    %%%
    %%% The script assumes that you have determined the z-offsets previously and
    %%% loads them z-offsets from a file
    %%%
    %%% Author: Jan Lipfert
    %%% Date: 2013-09-16
    
    clear all; clc; close all;
    
    %Datafile and aquisition frequency
    traces_file = 'DaldropData.txt';
    freq = 2800; %%% Acquisition frequency in Hz
    
    
    %Constants 
    kT = 4.1; %pN nm
    eta = 10E-10; %viscosity in pN s/nm^2
    Rbead = 515; % Bead radius in nm
  
    %%%--- Read in data ---
    data = load(traces_file);    
    Nbeads = 1;
    
    %%%--- Parse the bead data ---
    for i=1:Nbeads
        display(['Parsing data for bead ' num2str(i) ' of ' num2str(Nbeads)])
        bead(i).time = 1:length(data(:,1));
        bead(i).x = data(:,1)*1000;
        bead(i).y = data(:,2)*1000;
        bead(i).z = data(:,3)*1000; %nm
        
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Top level section 3:  FX3_ANALYZE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1
    
    %%% This is the third section in a series of routines to analyze
    %%% force-extension data taken with the new multi bead code for MT
    %%% (J. P. Cnossen, D. Dulin and N. H. Dekker Rev Sci Instrum 85, (2014).)
    %%%
    %%% This section assumes that you have read the (z-offset corrected) raw
    %%% data traces
    %%%
    %%% Author: Jan Lipfert
    %%% Date: 2013-09-16, updated 2014-07-17
    
    
    clc;
    
    %%% ------------------------------------------------------- %%%
    %%% Look over the beads and determine the force for each plateau
    %%% ------------------------------------------------------- %%%
    if 1
        for i=1:Nbeads
            
            [Ext, Fx_real, Fy_real, PSDfit, PSDforce, fcorner, MLfitx, MLforcex, Rfitx, MLfity, MLforcey, Rfity]=...
                analyze_one_trace2(bead(i).time,...
                bead(i).x,... 
                bead(i).y,...
                bead(i).z,...
                freq, Rbead, kT, eta);

            display(['Finished working on bead # ' num2str(i)])
            bead(i).ext = Ext;
            bead(i).Fx_real = Fx_real;
            bead(i).Fy_real = Fy_real;
            PSDfit = PSDfit;
            bead(i).PSDforce = PSDforce;
            bead(i).fcorner = fcorner;
            PSD2fitx = MLfitx;
            bead(i).PSD2force = MLforcex;
            bead(i).Rfitx = Rfitx;
            PSD3fity = MLfity;
            bead(i).PSD3force = MLforcey;
            bead(i).Rfity = Rfity;
            
            
            %calculate fc for y direction, just to check
            Cpar = (1-9/16*(1+(Ext)/Rfity)^(-1)+1/8*(1+(Ext)/Rfity)^(-3)-45/256*(1+(Ext)/Rfity)^(-4)-1/16*(1+(Ext)/Rfity)^(-5))^(-1); %Daldrop eq(S10)
            alphaY = 6*pi*eta*Rfity*Cpar; %Daldrop eq(11)
            kappa = MLforcey/(Ext);

            fc = kappa/(2*pi*alphaY);
            
            %calculate fmin and fplus for x direction, also to check
            Cpar = (1-9/16*(1+(Ext)/Rfitx)^(-1)+1/8*(1+(Ext)/Rfitx)^(-3)-45/256*(1+(Ext)/Rfitx)^(-4)-1/16*(1+(Ext)/Rfitx)^(-5))^(-1); %Daldrop eq(S10)
            Crot = 1 + 5/16*(1+(Ext)/Rfitx)^(-3);
            alphaX = 6*pi*eta*Rfitx*Cpar + 8*pi*eta*Rfitx*Crot/(1+Ext/Rfitx)^2; %Daldrop eq(11)
            alphaPhi = 8*pi*eta*Rfitx^3*Crot; %Daldrop eq(11)
                        
            fPlus = (MLforcex/(Ext)*((Ext+Rfitx)*Rfitx/(2*alphaPhi) + 1/(2*alphaX) + 1/2*(((Ext+Rfitx)*Rfitx/alphaPhi + 1/alphaX)^2-4*Ext*Rfitx/(alphaX*alphaPhi))^(1/2)))/(2*pi);
            fMin = (MLforcex/(Ext)*((Ext+Rfitx)*Rfitx/(2*alphaPhi) + 1/(2*alphaX) - 1/2*(((Ext+Rfitx)*Rfitx/alphaPhi + 1/alphaX)^2-4*Ext*Rfitx/(alphaX*alphaPhi))^(1/2)))/(2*pi);
            display(['Zhongbo method: F = ' num2str(bead(1).PSDforce) 'pN, R = ' num2str(Rbead) 'nm'])
            display(['Daldrop method y-direction: F = ' num2str(bead(1).PSD3force) 'pN, R = ' num2str(bead(1).Rfity) 'nm'])
            display(['Daldrop method x-direction: F = ' num2str(bead(1).PSD2force) 'pN, R = ' num2str(bead(1).Rfitx) 'nm'])
            display(['Labview implementation of Daldrop y-direction: F = 4.9115pN, R = 512nm'])
            display(['Labview implementation of Daldrop x-direction: F = 4.9223pN, R = 514nm'])
        end
    end
end