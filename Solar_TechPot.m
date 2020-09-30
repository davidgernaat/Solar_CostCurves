%% Technical potential CSP PV and PVres
KW=1;

SM = 2;               % Solar Multiple
CSPm2pMW=6000;        % amount of m2 needed for reference 1MW plant Trieb et al 2009
CSPm2pkW=6;           % amount of m2 needed for reference 1kW plant Trieb et al 2009

%% Temperature efficiency effect CSP
CSP_SolarElecEff 	= 0.12;	%Typical conversion efficiency for dry-cooled parabolic trough plant = 12% - Trieb et al 2009, Turchi et al 2010. Below we decompose this number to include climate impacts

% The efficiency of the CSP system depends on the optical efficiency of the mirrors and receiver, and the thermal efficiency of the receiver system and
% the heat transfer fluid. The optical efficiency is the percentage which determines what part of the incoming solar radiation which is absorbed by the
% receiver tube while the thermal efficiency is the percentage of that radiation transferred from the receiver to the heat transfer fluid and finally to the
% power block.

% Formule comes from: Crook et al 2011 Climate impacts PV and CSP
k0 = 0.762; % 
k1 = 0.2125; % W m-2 Celcius
Tref = 115; % dCelcius

CSP_ThermalEff = 0.3; %Rankine cycle losses
CSP_SystemAvail = 0.9; %System availability
for m=1:13
    CSP_EffTot{m} = zeros(nr,nc);
    for r=1:nr
        for c=1:nc
            if Avg_DNI_NASA_adj_Insol{m}(r,c)==0; CSP_ReceiverEff{m}(r,c)=0; continue; end
            CSP_ReceiverEff{m}(r,c) = max(0, k0 - ( (k1 * (Tref - tasc{m}(r,c))) / Avg_DNI_NASA_adj_Insol{m}(r,c)));
            CSP_EffTot{m}(r,c) = CSP_ReceiverEff{m}(r,c) * CSP_ThermalEff * CSP_SystemAvail;
        end
    end
    CSP_EffTot{m}(find(isinf(CSP_EffTot{m}(:))))=0; %remove infs
end

% figure(1);clf;
% ax1 = subplot(1,2,1); imagesc(PRcc{13});axis image; colorbar; title('Temperature')
% ax2 = subplot(1,2,2); imagesc(PRcc{14});axis image; colorbar; title('PR')
% linkaxes([ax1 ax2])

% figure(1);clf;imagesc(tasc{13});axis image
% figure(2);clf;imagesc(Avg_DNI_NASA_adj_Insol{13});axis image
% figure(3);clf;imagesc(CSP_EffTot{13});axis image; colorbar

%% DNI classes This is the old code from Alex

% CSP_LoadFac = zeros(nr,nc);
% 
% Avg_DNI_NASA_adj_ann = Avg_DNI_NASA_adj{13} * 365; %kWh / m2 / y
% 
% for r=1:nr
%     for c=1:nc
%         if Avg_DNI_NASA_adj_ann(r,c) < 1095
%             CSP_LoadFac(r,c) = 0;
%         elseif Avg_DNI_NASA_adj_ann(r,c) > 1095 && Avg_DNI_NASA_adj_ann(r,c) < 2800
%             
%             CSP_LoadFac(r,c) = 0.00022041 * Avg_DNI_NASA_adj_ann(r,c) -0.039;
%             
%         elseif Avg_DNI_NASA_adj_ann(r,c) > 2800
%             CSP_LoadFac(r,c) = 0.6;
%         end
%     end
% end
% figure(2);clf;imagesc(CSP_LoadFac);axis image; colorbar
% 
% mean(CSP_Loadfac(Cind{166}))

%% CSP loadfactor This is new and similar to the way PV load factor is calculated but multiplied with solar multiple
% The idea with solar multiple is that if you have 2x the solar field (mirrors)
% required for the turbine, you can store that thermal energy and use that
% energy in the evening, thus increasing your loadfactor with 2x.

for m=13
    for r=1:nr
        for c=1:nc
            if Avg_DNI_NASA_adj_Insol{m}(r,c)==0; CSP_LoadFac(r,c)=0; continue; end;
            if CSP_EffTot{m}(r,c)==0; CSP_LoadFac(r,c)=0; continue; end;
                
            CSPConvWperm2{m}(r,c) = 1/((CSP_EffTot{m}(r,c))*1000);
            CSPWPM2{m}(r,c) = 1/CSPConvWperm2{m}(r,c); %Amount of watt per m2 depending on efficiencies
            CSP_LoadFac(r,c) = ((Avg_DNI_NASA_adj_Insol{m}(r,c) * CSP_EffTot{m}(r,c)) / CSPWPM2{m}(r,c)) * SM ; % Load factor per cell
            
        end
    end
end

% figure(1);clf;imagesc(CSP_LoadFac);axis image; colorbar
% mean(CSP_Loadfac2(Cind{166}))
            
%% CSP Technical potential per m2 per cell
% It is assumed, based on Trieb et al., that everything is in the efficiency, also the loadfactor.

for r=1:nr
    for c=1:nc
        if (Avg_DNI_NASA_adj{13}(r,c)*365) < 1095 %minimum amount of irradiation Trieb et al
            CSP_TechPot_Cell_LF(r,c)=0;
        else
            CSP_TechPot_Cell_LF(r,c) = Avg_DNI_NASA_adj{13}(r,c) * 365 * CSP_EffTot{13}(r,c); % kWh / m2 / y |  used in cost-curves
        end        
    end
end

% figure(1);clf;imagesc(CSP_TechPot_Cell_LF);axis image;colorbar
% figure(2);clf;imagesc(Avg_DNI_NASA_adj_ann);axis image;colorbar

%% CSP Tech Regional and Global Potentials

CSP_LandUseFactor 	= 0.37;	%Mid-range value from Trieb et al 2009 = 37% of the land is covered by solar reflectros

for r=1:nr
    for c=1:nc
        AnnTechPotCellCSP(r,c) = CSP_GeoPot_Cell{13}(r,c) ...
            * CSP_LandUseFactor * CSP_EffTot{m}(r,c) * CSP_LoadFac(r,c); % kWh / cell / y
    end
end

% figure(1);clf;imagesc(AnnTechPotCellCSP);colorbar

for i=1:26
    RegTechPotCSP(i) = sum(AnnTechPotCellCSP(IRind{i})); % kWh / region / y
end
RegTechPotCSP(27) =sum(RegTechPotCSP(1:26));

for i=1:numel(ISOGDP(:,1))
    CTechPotCSP(i) = sum(AnnTechPotCellCSP(Cind{i})); % kWh / country / y
end
dt = sum(CTechPotCSP(1:end));
CTechPotCSP(end+1) = dt;

%% PV

%% Temperature map for climate effect on PV efficiency
% Formule comes from Jerez et al 2014
c1 = 4.3; %dCelcius
c2 = 0.943; %
c3 = 0.028; %dCelcius m2 W-1
c4 = -1.528; %dCelcius s m-1
for m=1:12
    for r=1:nr
        for c=1:nc
            T{m}(r,c) = c1 + c2*tasc{m}(r,c) + c3*Gl_Horiz_NASA_Insol{m}(r,c) + c4*wsc{m}(r,c); %Temperature map (C)
            PRcc{m}(r,c) = min(1,1 + -0.005*(T{m}(r,c) - 25)); %Change in efficiency through climate variables
        end
    end
end

% Yearly average efficiency based on months. Otherwise no efficiency loss
% in Europe for example. Better still would be yearly average based on days.
for r=1:nr
    for c=1:nc
        clear dt
        for m=1:12
            dt(m) = PRcc{m}(r,c);
        end
        PRcc{13}(r,c) = mean(dt);
    end
end

% figure(1);clf;
% for m=1:13
%     ax{m} = subplot(4,4,m); imagesc(PRcc{m});axis image; colorbar; title(sprintf('Month %d',m))
% end
% linkaxes([ax{:}])

% figure(1);clf;
% ax1 = subplot(1,2,1); imagesc(PRcc{13});axis image; colorbar; title('Temperature')
% ax2 = subplot(1,2,2); imagesc(PRcc{14});axis image; colorbar; title('PR')
% linkaxes([ax1 ax2])

% PV efficiency effect
PV_SolarElecEff   = 0.17; % IRENA 2018 Renewble power costs and Fraunhofer. 2018. Photovoltaics report, Fraunhofer Institute for Solar Energy Systems
Performance_ratio = 0.85; % Fraunhofer. 2018. Photovoltaics report, Fraunhofer Institute for Solar Energy Systems
PVeff_nocc        = PV_SolarElecEff * Performance_ratio;

for m=13
    for r=1:nr
        for c=1:nc
            PVeff{m}(r,c) = PV_SolarElecEff * Performance_ratio * PRcc{m}(r,c);
            
            ConvWperm2{m}(r,c) = 1/((PVeff{m}(r,c))*1000);
            WPM2{m}(r,c) = 1/ConvWperm2{m}(r,c); %Amount of watt-peak per m2 depending on efficiencies

        end
    end
end

% figure(1);clf;imagesc(PVeff{13});axis image; colorbar
% figure(2);clf;imagesc(WPM2{13});axis image; colorbar

%% Loadfactors

for m=13
    for r=1:nr
        for c=1:nc
            PV_Loadfac_Cell{m}(r,c) = (Gl_Horiz_NASA_Insol{m}(r,c) * PVeff{m}(r,c)) / WPM2{m}(r,c) ; % Load factor per cell
            
            PV_TechPot_Cell_LF{m}(r,c) = WPM2{m}(r,c)*1e-3 * 8760 * PV_Loadfac_Cell{m}(r,c); % kWh / m2 / y
            
            PV_TechPot_Cell_Eff{m}(r,c) = 365 * Gl_Horiz_NASA{m}(r,c) * PVeff{m}(r,c); % kWh / m2 / y
        end
    end
end
% figure(1);clf;imagesc(PV_Loadfac_Cell{13});axis image;colorbar

%% PV Tech Regional and Global Potentials
PV_LandUseFactor 	= 0.47;	%PV Average Land Use coverage by panels (Packing Factor reported in Ong et al 2013)

for r=1:nr
    for c=1:nc
        AnnTechPotCellPV(r,c) = PV_GeoPot_Cell{13}(r,c) ...
            * PVeff{13}(r,c) * PV_LandUseFactor; % kWh / cell / y
    end
end

for i=1:26
    RegTechPotPV(i) = sum(AnnTechPotCellPV(IRind{i})); % kWh / region / y
end
RegTechPotPV(27) =sum(RegTechPotPV(1:26));

for i=1:numel(ISOGDP(:,1))
    CTechPotPV(i) = sum(AnnTechPotCellPV(Cind{i})); % kWh / country / y
end
dt = sum(CTechPotPV(1:end));
CTechPotPV(end+1) = dt;


%% PVres
PVres_Loadfac_Cell= zeros(nr,nc);
PVres_Loadfac_Cell_LF= zeros(nr,nc);
PVres_Loadfac_Cell_Eff= zeros(nr,nc);

for m=13;
    for r=1:nr
        for c=1:nc
            PVres_Loadfac_Cell(r,c) = (Gl_Horiz_NASA_Insol{13}(r,c) * PVeff{m}(r,c)) / WPM2{m}(r,c) ; % Load factor per cell
            
            PVres_TechPot_Cell_LF(r,c) = WPM2{m}(r,c)*1e-3 * 8760 * PVres_Loadfac_Cell(r,c); % kWh / m2 / y | used in cost-curve | efficiency is already in WPM2
            
            PVres_TechPot_Cell_Eff(r,c) = 365 * Gl_Horiz_NASA{m}(r,c) * PVeff{m}(r,c); % kWh / m2 / y
        end
    end
end

%% PVres Tech Regional and Global Potentials

AnnTechPotCellPVres = zeros(nr,nc);
for m=13;
    for r=1:nr
        for c=1:nc
            AnnTechPotCellPVres(r,c) = PVres_GeoPot{m}(r,c) ... %Total
                * PVeff{m}(r,c);
            
            AnnTechPotCellPVresu(r,c) = PVres_GeoPotu{m}(r,c) ... %Urban
                * PVeff{m}(r,c);
            
            AnnTechPotCellPVresr(r,c) = PVres_GeoPotr{m}(r,c) ... %Rural
                * PVeff{m}(r,c);
            
        end
    end
end

for i=1:26
    RegTechPotPVres(i) = sum(AnnTechPotCellPVres(IRind{i})); % kWh / region / y
    RegTechPotPVresu(i) = sum(AnnTechPotCellPVresu(IRind{i})); % kWh / region / y Urban
    RegTechPotPVresr(i) = sum(AnnTechPotCellPVresr(IRind{i})); % kWh / region / y Rural
    
end
RegTechPotPVres(27) =sum(RegTechPotPVres(1:26));
RegTechPotPVresu(27) =sum(RegTechPotPVresu(1:26));
RegTechPotPVresr(27) =sum(RegTechPotPVresr(1:26));

% fprintf('Global: %.0f TWh\n',(RegTechPotPVres_nocc(27))*1e-9)
% fprintf('Global: %.0f TWh\n',(RegTechPotPVresu_nocc(27))*1e-9)
% fprintf('USA: %.0f TWh\n',(RegTechPotPVres_nocc(2))*1e-9)
% fprintf('India: %.0f TWh\n',(RegTechPotPVres_nocc(18))*1e-9)
% fprintf('EU: %.0f TWh\n',(RegTechPotPVres_nocc(11)+RegTechPotPVres_nocc(12))*1e-9)
% fprintf('Brazil: %.0f TWh\n',(RegTechPotPVres_nocc(5))*1e-9*0.15)
% fprintf('USA without roofuse factor: %.2f TWh\n',(RegTechPotPVres_nocc(2)*(1/PVres_RoofUseFactor))*1e-9)

for i=1:numel(ISOGDP(:,1))
    CTechPotPVres(i) = sum(AnnTechPotCellPVres(Cind{i})); % kWh / country / y
    CTechPotPVresu(i) = sum(AnnTechPotCellPVresu(Cind{i})); % kWh / country / y
    CTechPotPVresr(i) = sum(AnnTechPotCellPVresr(Cind{i})); % kWh / country / y
    
end
dt = sum(CTechPotPVres(1:end));
CTechPotPVres(end+1) = dt;
dt = sum(CTechPotPVresu(1:end));
CTechPotPVresu(end+1) = dt;
dt = sum(CTechPotPVresr(1:end));
CTechPotPVresr(end+1) = dt;

% fprintf('Switserland %.0f TWh\n',CTechPotPVres_nocc(171)*1e-9);
% fprintf('Spain %.0f TWh\n',CTechPotPVres_nocc(166)*1e-9);

%% Nederland
% a=PV_TechPot_Cell_LF{13};
% a(find(GISO(:)==528)) %Nederland
% figure(1);clf;imagesc(PV_TechPot_Cell_LF{13});colorbar
% 
% b=AnnTechPotCellPVres./RAmap;
% 
% c=AnnTechPotCellPVres(find(GISO(:)==528));
% d=RAmap(find(GISO(:)==528));
% e=(c./d)/0.32
% % f=PVeff{m}(find(GISO(:)==528))
% 
% figure(2);clf;imagesc(PVres_Loadfac_Cell);colorbar

%% Specific area
% win=0;
% a = IRegion==11;
% cols = find(sum(a)>0);
% rows = find(sum(a')>0);
% fc = cols(1) - win;
% lc = cols(end) + win;
% fr = rows(1) - win;
% lr = rows(end) + win;
% 
% b=PR{13};
% b(find(~a))=0.9;
% 
% bw = b(fr:lr,fc:lc);
% 
% figure(1);clf;imagesc(bw);colorbar; axis image
% % figure(2);clf;hist(PR{13}(:))

%% PV potential Germany
% figure(1);clf;imagesc(GISO);colorbar; axis image
% a = GISO==276;
% b=AnnTechPotCellPV;
% b(find(~a))=0;
% figure(1);clf;imagesc(b);colorbar; axis image
% 
% fprintf('Solar potential: %0.2f TWh\n',sum(b(:))*1e-9)
% fprintf('Current generation Germany: %0.2f TWh\n',46.164)
