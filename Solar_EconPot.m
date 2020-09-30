%%Calculation the cost per cell

Interest = 0.1;         % 10
EconLifeTime = 20;      % 20 years for all technologies
OnMShare = 0.03;        % 3%	Reference: NEEDS 2008, p.29
AnnualHours	= 8760;     % Number of hours in a year
NewTransmissionCost_perKM = 1843000; % Ontario case study. Veileder is 2x cheaper
% NewTransmissionCost_perKM = 2390000; % Cost per km of building a new 230 kV double circuit transmission line (Source: Mason et al 2012)
AnnuityFactor = Interest/(1-(1+Interest)^-EconLifeTime);

PV_SpecInvest = 1288;    % $2010/kW See PVres paper Europe %Fraunhofer 2018 rooftop pv:1400 EUR/kW; larger systems 1140 €/kWp=1288 $/kWh
PVres_SpecInvest = 280;  % $2010/kW See PVres paper
CSP_SpecInvest = 7000;   % $2010/kW IRENA Cost analysis for CSP plant 6hrs storage p1 in 2010

CSP_OnMCost		= 70;       % $/KW-yr for PT reference plant, (NREL, Turchi et al 2010)
PV_OnMCost		= 16;		% $/KW-yr See PVres paper Europe before:for PV reference plant, (estimeted at 2% of Installed costs by Desideri and Campana 2014)
PVres_OnMCost	= 16;		% $/KW-yr See PVres paper Europe before:for PVres plant, (NREL 2013, Distributed Generation Renewable Energy Estimate of Costs, accessed 29.03.2016)

DistanceToLoad2 = DistanceToLoad;
DistanceToLoad2(DistanceToLoad>5000)=5000;

%% CSP Capacity per cell
for r=1:nr
    for c=1:nc
        CSPCap(r,c) = (Area(r,c) * 1e6 * SuitabilityFactor(r,c) * BioMult(r,c) * BuildupInv(r,c) * Slope(r,c)) *0.37 ... % Remove all unsuitable area
            / CSPm2pMW; % Divided by the amount of m2 per 1MW CSP plant
        % CSP plants / cell or MW / cell
    end
end

%% CSP Distance to Load powerline cost
for r=1:nr
    for c=1:nc
        CSPDisLoadCost(r,c) = Powerline_allocator(DistanceToLoad2(r,c),CSPCap(r,c),0) ; % $2010 inputs Distance km, Capcity MW, 0 = no multiple powerlines
        CSPDisLoadCostKW(r,c) = CSPDisLoadCost(r,c) / (CSPCap(r,c)*1e3);
        % Remove infs
        if isinf(CSPDisLoadCostKW(r,c))==1; CSPDisLoadCostKW(r,c)=NaN; end
    end
end

%% CSP Cost

CSP_TotalAnnualCost = KW * AnnuityFactor * CSP_SpecInvest + CSP_OnMCost;	% total Annual Cost for 1 kW

% Adding powerline cost
for r=1:nr
    for c=1:nc
        
        CSPTransCostpKW(r,c) = (DistanceToLoad2(r,c) * NewTransmissionCost_perKM) / (CSPCap(r,c)*1e3); %$ / kW
        
        % Remove infs
        if isinf(CSPTransCostpKW(r,c))==1; CSPTransCostpKW(r,c)=NaN; end
        
        CSP_TotalAnnualCost_PL(r,c) = CSP_TotalAnnualCost + AnnuityFactor * CSPTransCostpKW(r,c); % $
        
        CSP_TotalAnnualCost_PLDL(r,c) = CSP_TotalAnnualCost + AnnuityFactor * CSPDisLoadCostKW(r,c); % $
    end
end

%% CSP COE

for r=1:nr
    for c=1:nc
        CSP_COE_LF(r,c) = (CSP_TotalAnnualCost/(CSPm2pkW*SM)) / CSP_TechPot_Cell_LF(r,c) ; % $ / kWh | cost is per kW so needs to be converted to m2
        
        CSP_COE_PL(r,c) = (CSP_TotalAnnualCost_PL(r,c)/(CSPm2pkW*SM)) / CSP_TechPot_Cell_LF(r,c); % $ / kWh | used in the cost-curves
        CSP_COE_PLDL(r,c) = (CSP_TotalAnnualCost_PLDL(r,c)/(CSPm2pkW*SM)) / CSP_TechPot_Cell_LF(r,c); % $ / kWh
    end
end
CSP_COE_LF(find(isinf(CSP_COE_LF(:))))=NaN; %remove infs
CSP_COE_PL(find(isinf(CSP_COE_PL(:))))=NaN; %remove infs
CSP_COE_PLDL(find(isinf(CSP_COE_PLDL(:))))=NaN; %remove infs

% CSP_COE_LF(CSP_COE_LF(:)>1)=1;
% CSP_COE_PL(CSP_COE_PL>0.5)=0.5;
% CSP_COE_PLDL(CSP_COE_PLDL>0.5)=0.5;
% figure(1);clf;
% ax1 = subplot(2,1,1);
% imagesc(CSP_TechPot_Cell_LF./8760); axis image; colorbar;title('CSP')
% ax2 = subplot(2,1,2);
% imagesc(AnnTechPotCellCSP); axis image; colorbar;
% linkaxes([ax1 ax2],'xy')

% figure(1);clf;imagesc(CSP_TechPot_Cell_LF);axis image;colorbar
% a=CSP_COE_PL;
% a(a(:)>1)=1;
% figure(1);clf;imagesc(CSP_EffTot{13});axis image;colorbar
% figure(2);clf;imagesc(a);axis image;colorbar

%% Remove all the cells that have zero potential as defined in GeoPot

CSP_LoadFac_prep = CSP_LoadFac;
for r=1:nr
    for c=1:nc
        if CSP_GeoPot_Cell{13}(r,c)==0;
            CSP_COE_LF(r,c)=inf;
            CSP_COE_PL(r,c)=inf;
            CSP_COE_PLDL(r,c) = inf;
            CSP_LoadFac_prep(r,c)=inf;
        end
    end
end

if setEconCeil==1
    CSP_COE_LF(CSP_COE_LF>CostCeil)=NaN;
    CSP_COE_PL(CSP_COE_PL>CostCeil)=NaN;
    CSP_COE_PLDL(CSP_COE_PLDL>CostCeil)=NaN;
    if setCSPPL==0
        CSP_LoadFac_prep(CSP_COE_LF>CostCeil)=NaN;
    elseif setCSPPL==1
        if multifactor==0
            CSP_LoadFac_prep(CSP_COE_PL>CostCeil)=NaN;
        elseif multifactor==1
            CSP_LoadFac_prep(CSP_COE_PLDL>CostCeil)=NaN;
        end
    end
end

% figure(1);clf;imagesc(CSP_COE_Eff);axis image; colorbar
% figure(2);clf;imagesc(CSP_COE_LF);axis image; colorbar

%% CSP Cost curve prep per IMAGE region
clear SI
%Find cell per Image reg
for i=1:26
    
    if setCSPPL==0
        CSP_COE_IR{i} = CSP_COE_LF(IRind{i});
    elseif setCSPPL==1
        if multifactor==0
            CSP_COE_IR{i} = CSP_COE_PL(IRind{i});
        elseif multifactor==1
            CSP_COE_IR{i} = CSP_COE_PLDL(IRind{i});
        end
    end
    
    CSP_Pot_IR{i} = AnnTechPotCellCSP(IRind{i});
    
    CSP_LF_IR{i} = CSP_LoadFac_prep(IRind{i});
    
end

for i=1:(nr*nc)
    if setCSPPL==0
        CSP_COE_IR{27}(i) = CSP_COE_LF(i);
    elseif setCSPPL==1
        if multifactor==0
            CSP_COE_IR{27}(i) = CSP_COE_PL(i);
        elseif multifactor==1
            CSP_COE_IR{27}(i) = CSP_COE_PLDL(i);
        end
    end
    
    CSP_Pot_IR{27}(i) = AnnTechPotCellCSP(i);
    CSP_LF_IR{27}(i) = CSP_LoadFac_prep(i);
end

for i=1:26
    [CSP_COE_IR_S{i}, SI{i}] = sort(CSP_COE_IR{i});
    CSP_Pot_IR_S{i} = CSP_Pot_IR{i}(SI{i});
    CSP_LF_IR_S{i} = CSP_LF_IR{i}(SI{i});
end

[CSP_COE_IR_S{27}, SI{27}] = sort(CSP_COE_IR{27});
CSP_Pot_IR_S{27} = CSP_Pot_IR{27}(SI{27});
CSP_LF_IR_S{27} = CSP_LF_IR{27}(SI{27});

for i=1:27
    CSP_COE_IR_S{i} = CSP_COE_IR_S{i}(~isnan(CSP_COE_IR_S{i}));
    CSP_Pot_IR_S{i} = CSP_Pot_IR_S{i}(~isnan(CSP_COE_IR_S{i}));
    CSP_LF_IR_S{i} = CSP_LF_IR_S{i}(~isnan(CSP_LF_IR_S{i}));
    
    CSP_COE_IR_S{i} = CSP_COE_IR_S{i}(find(~isinf(CSP_COE_IR_S{i})));
    CSP_Pot_IR_S{i} = CSP_Pot_IR_S{i}(find(~isinf(CSP_COE_IR_S{i})));
    CSP_LF_IR_S{i} = CSP_LF_IR_S{i}(find(~isinf(CSP_LF_IR_S{i})));
    
end

for i=1:27
    RegEconPotCSP(i) = sum(CSP_Pot_IR_S{i}); % kWh / region / y
end;

%% CSP Cost curve prep per country
clear SI

ncw=(numel(ISOGDP(:,1))+1); %Number of countries plus world

%Find cell per country
for i=1:(ncw-1)
    if setCSPPL==0
        CSP_COE_C{i} = CSP_COE_LF(Cind{i});
    elseif setCSPPL==1
        if multifactor==0
            CSP_COE_C{i} = CSP_COE_PL(Cind{i});
        elseif multifactor==1
            CSP_COE_C{i} = CSP_COE_PLDL(Cind{i});
        end
    end
    
    CSP_Pot_C{i} = AnnTechPotCellCSP(Cind{i});
    
    CSP_LF_C{i} = CSP_LoadFac_prep(Cind{i});
    
end

for i=1:(nr*nc)
    if setCSPPL==0 %No powerline cost included
        CSP_COE_C{ncw}(i) = CSP_COE_LF(i);
    elseif setCSPPL==1 %Powerline cost included
        if multifactor==0 
            CSP_COE_C{ncw}(i) = CSP_COE_PL(i);
        elseif multifactor==1 %Multiple powerline cost included
            CSP_COE_C{ncw}(i) = CSP_COE_PLDL(i);
        end
    end
    
    CSP_Pot_C{ncw}(i) = AnnTechPotCellCSP(i);
    CSP_LF_C{ncw}(i) = CSP_LoadFac_prep(i);
end

for i=1:(ncw-1)
    [CSP_COE_C_S{i}, SI{i}] = sort(CSP_COE_C{i});
    CSP_Pot_C_S{i} = CSP_Pot_C{i}(SI{i});
    CSP_LF_C_S{i} = CSP_LF_C{i}(SI{i});
end

[CSP_COE_C_S{ncw}, SI{ncw}] = sort(CSP_COE_C{ncw});
CSP_Pot_C_S{ncw} = CSP_Pot_C{ncw}(SI{ncw});
CSP_LF_C_S{ncw} = CSP_LF_C{ncw}(SI{ncw});

for i=1:ncw
    CSP_COE_C_S{i} = CSP_COE_C_S{i}(~isnan(CSP_COE_C_S{i}));
    CSP_Pot_C_S{i} = CSP_Pot_C_S{i}(~isnan(CSP_COE_C_S{i}));
    CSP_LF_C_S{i} = CSP_LF_C_S{i}(~isnan(CSP_LF_C_S{i}));
    
    CSP_COE_C_S{i} = CSP_COE_C_S{i}(find(~isinf(CSP_COE_C_S{i})));
    CSP_Pot_C_S{i} = CSP_Pot_C_S{i}(find(~isinf(CSP_COE_C_S{i})));
    CSP_LF_C_S{i} = CSP_LF_C_S{i}(find(~isinf(CSP_LF_C_S{i})));
    
end

for i=1:ncw
    CEconPotCSP(i) = sum(CSP_Pot_C_S{i}); % kWh / country / y
end

%% PV Capacity per cell
m=13;
for r=1:nr
    for c=1:nc
        PVCap(r,c) = WPM2{m}(r,c) ... % W / m2 PV paneel
            * (Area(r,c) * 1e6) ... % Maal het aantal m2 per cell
            * SuitabilityFactor(r,c) * BioMult(r,c) * BuildupInv(r,c) * Slope(r,c) ... % remove all unsuitable area
            * PV_LandUseFactor * 1e-6; % Include landuse efficiency and convert to MW
        % MW / cell
    end
end

%% PV DisCost
for r=1:nr
    for c=1:nc
        PVDisLoadCost(r,c) = Powerline_allocator(DistanceToLoad2(r,c),PVCap(r,c),0) ; % $2010 inputs Distance km, Capcity MW, 0 = no multiple powerlines
        PVDisLoadCostKW(r,c) = PVDisLoadCost(r,c) / (PVCap(r,c)*1e3);
    end
end

%% PV Cost

PV_TotalAnnualCost = AnnuityFactor * PV_SpecInvest * WPM2{m}(r,c)*1e-3 + (PV_OnMCost * WPM2{m}(r,c)*1e-3);	% Total Annual Cost for 1 KW

% Adding powerline cost
for r=1:nr
    for c=1:nc
        
        PVTransCostpKW(r,c) = (DistanceToLoad2(r,c) * NewTransmissionCost_perKM) / (Area(r,c)*1e6*WPM2{m}(r,c)*1e-3*PV_LandUseFactor*0.1); %$ / kW %Last factor is a supposedly suitability factor
        PV_TotalAnnualCost_PL(r,c) = PV_TotalAnnualCost + PVTransCostpKW(r,c)*AnnuityFactor; % $
        PV_TotalAnnualCost_PLDL(r,c) = PV_TotalAnnualCost + PVDisLoadCostKW(r,c)*AnnuityFactor; % $
        
        PVTransCostpKWh(r,c) = (DistanceToLoad2(r,c) * NewTransmissionCost_perKM) / AnnTechPotCellPV(r,c);
        PVTransCostpKWhDL(r,c) = PVDisLoadCost(r,c) / AnnTechPotCellPV(r,c);
        
    end
end

%% PV COE
for r=1:nr
    for c=1:nc
        PV_COE_Eff(r,c) = PV_TotalAnnualCost / PV_TechPot_Cell_Eff{13}(r,c); % $ / kWh
        PV_COE_LF(r,c) = PV_TotalAnnualCost / PV_TechPot_Cell_LF{13}(r,c); % $ / kWh
        PV_COE_PL(r,c) = PV_TotalAnnualCost_PL(r,c) / PV_TechPot_Cell_LF{13}(r,c); % $ / kWh
        PV_COE_PLDL(r,c) = (PV_TotalAnnualCost_PLDL(r,c) / PV_TechPot_Cell_LF{13}(r,c)); % $ / kWh
        
        PV_COE_PLh(r,c) = (PV_TotalAnnualCost / PV_TechPot_Cell_LF{13}(r,c)) + PVTransCostpKWh(r,c); % $ / kWh | used in cost-curve
        PV_COE_PLhDL(r,c) = (PV_TotalAnnualCost / PV_TechPot_Cell_LF{13}(r,c)) + PVTransCostpKWhDL(r,c); % $ / kWh
        
    end
end

% PV_COE_PLh(PV_COE_PLh>0.5)=0.5;
% PV_COE_PLhDL(PV_COE_PLhDL>0.5)=0.5;
% figure(2);clf;
% ax1 = subplot(2,1,1);
% imagesc(PV_COE_PLh); axis image; colorbar;title('PV')
% ax2 = subplot(2,1,2);
% imagesc(PV_COE_PLhDL); axis image; colorbar;
% linkaxes([ax1 ax2],'xy')

% PVTransCostpKWh(isinf(PVTransCostpKWh))=0.5;
% PVTransCostpKWh(PVTransCostpKWh>1)=1;
% figure(1);clf;imagesc(PVTransCostpKWh);axis image
% 
% figure(3);clf;imagesc(DistanceToLoad2);axisi image; colorbar

%% Remove all the cells that have zero potential as defined in GeoPot
PV_Loadfac_Cell_prep = PV_Loadfac_Cell{13};

for r=1:nr
    for c=1:nc
        if PV_GeoPot_Cell{13}(r,c)==0;
            PV_COE_LF(r,c)=NaN;
            PV_COE_Eff(r,c)=NaN;
            PV_COE_PL(r,c)=NaN;
            PV_COE_PLh(r,c)=NaN;
            PV_COE_PLDL(r,c)=NaN;
            
            PV_Loadfac_Cell_prep(r,c)=NaN;
            
        end
    end
end

if setEconCeil==1
    PV_COE_LF(PV_COE_LF>CostCeil)=NaN;
    PV_COE_Eff(PV_COE_Eff>CostCeil)=NaN;
    PV_COE_PL(PV_COE_PL>CostCeil)=NaN;
    PV_COE_PLh(PV_COE_PLh>CostCeil)=NaN;
    PV_COE_PLDL(PV_COE_PLDL>CostCeil)=NaN;
    if setPVPL==0
        PV_Loadfac_Cell_prep(PV_COE_LF>CostCeil)=NaN;
    elseif setPVPL==1
        if multifactor==0
            PV_Loadfac_Cell_prep(PV_COE_PLh>CostCeil)=NaN;
        elseif multifactor==1
            PV_Loadfac_Cell_prep(PV_COE_PLDL>CostCeil)=NaN;
        end
    end
end

% figure(3);clf;imagesc(PV_COE_LF);axis image; colorbar
% a=PV_COE_PLh;
% a(a>0.5)=0.5;
% figure(4);clf;imagesc(a);axis image; colorbar

%% PV Cost curve IMAGE regions
clear SI
%Find cell per Image reg
for i=1:26
    
    if setPVPL==0
        PV_COE_IR{i} = PV_COE_LF(IRind{i});
    elseif setPVPL==1
        if multifactor==0
            PV_COE_IR{i} = PV_COE_PLh(IRind{i});
        elseif multifactor==1
            PV_COE_IR{i} = PV_COE_PLDL(IRind{i});
        end
    end
    
    PV_Pot_IR{i} = AnnTechPotCellPV(IRind{i});
    
    PV_LF_IR{i} = PV_Loadfac_Cell_prep(IRind{i});
end

for i=1:(nr*nc)
    if setPVPL==0
        PV_COE_IR{27}(i) = PV_COE_LF(i);
    elseif setPVPL==1
        if multifactor==0
            PV_COE_IR{27}(i) = PV_COE_PLh(i);
        elseif multifactor==1
            PV_COE_IR{27}(i) = PV_COE_PLDL(i);
        end
    end
    
    PV_Pot_IR{27}(i) = AnnTechPotCellPV(i);
    PV_LF_IR{27}(i)  = PV_Loadfac_Cell_prep(i);
end
for i=1:26
    [PV_COE_IR_S{i}, SI{i}] = sort(PV_COE_IR{i});
    PV_Pot_IR_S{i} = PV_Pot_IR{i}(SI{i});
    PV_LF_IR_S{i} = PV_LF_IR{i}(SI{i});
end

[PV_COE_IR_S{27}, SI{27}] = sort(PV_COE_IR{27});
PV_Pot_IR_S{27} = PV_Pot_IR{27}(SI{27});
PV_LF_IR_S{27} = PV_LF_IR{27}(SI{27});

for i=1:27
    PV_COE_IR_S{i} = PV_COE_IR_S{i}(~isnan(PV_COE_IR_S{i}));
    PV_Pot_IR_S{i} = PV_Pot_IR_S{i}(~isnan(PV_COE_IR_S{i}));
    PV_LF_IR_S{i}  = PV_LF_IR_S{i}(~isnan(PV_COE_IR_S{i}));
end

for i=1:27
    RegEconPotPV(i) = sum(PV_Pot_IR_S{i}); % kWh / region / y
end;

%% PV Cost curve Countries
clear SI

for i=1:(ncw-1)
    
    if setPVPL==0
        PV_COE_C{i} = PV_COE_LF(Cind{i});
    elseif setPVPL==1
        if multifactor==0
            PV_COE_C{i} = PV_COE_PLh(Cind{i});
        elseif multifactor==1
            PV_COE_C{i} = PV_COE_PLDL(Cind{i});
        end
    end
    
    PV_Pot_C{i} = AnnTechPotCellPV(Cind{i});
    
    PV_LF_C{i} = PV_Loadfac_Cell_prep(Cind{i});
end

for i=1:(nr*nc)
    if setPVPL==0
        PV_COE_C{ncw}(i) = PV_COE_LF(i);
    elseif setPVPL==1
        if multifactor==0
            PV_COE_C{ncw}(i) = PV_COE_PLh(i);
        elseif multifactor==1
            PV_COE_C{ncw}(i) = PV_COE_PLDL(i);
        end
    end
    
    PV_Pot_C{ncw}(i) = AnnTechPotCellPV(i);
    PV_LF_C{ncw}(i)  = PV_Loadfac_Cell_prep(i);
end

for i=1:(ncw-1)
    [PV_COE_C_S{i}, SI{i}] = sort(PV_COE_C{i});
    PV_Pot_C_S{i} = PV_Pot_C{i}(SI{i});
    PV_LF_C_S{i} = PV_LF_C{i}(SI{i});
end

[PV_COE_C_S{ncw}, SI{ncw}] = sort(PV_COE_C{ncw});
PV_Pot_C_S{ncw} = PV_Pot_C{ncw}(SI{ncw});
PV_LF_C_S{ncw} = PV_LF_C{ncw}(SI{ncw});

for i=1:ncw
    PV_COE_C_S{i} = PV_COE_C_S{i}(~isnan(PV_COE_C_S{i}));
    PV_Pot_C_S{i} = PV_Pot_C_S{i}(~isnan(PV_COE_C_S{i}));
    PV_LF_C_S{i}  = PV_LF_C_S{i}(~isnan(PV_COE_C_S{i}));
end

for i=1:ncw
    CEconPotPV(i) = sum(PV_Pot_C_S{i}); % kWh / country / y
end;

%% PVres Cost

if regcapcost==1
    
    %From WEO 2016
    PVres_SpecInvest = [3480
        3480
        2680
        2680
        2680
        2680
        2840
        2840
        2840
        2840
        1600
        1600
        3480
        3480
        3480
        3480
        3000
        1460
        2880
        1480
        1480
        1480
        2880
        3480
        1480
        2840
        ];
    
    PVres_OnMCost = [34
        34
        30
        30
        30
        30
        28
        28
        28
        28
        16
        16
        34
        34
        34
        34
        30
        14
        28
        14
        14
        14
        28
        34
        14
        28];
    
    for i=1:numel(PVres_SpecInvest)
        PVres_TotalAnnualCost(i) = AnnuityFactor * (PVres_SpecInvest(i) + PVres_OnMCost(i)) * WPM2{m}(r,c)*1e-3 + (PVres_OnMCost(i) * WPM2{m}(r,c)*1e-3);	% Total Annual Cost for 1 KW
    end
    PVres_TotalAnnualCost_mean = mean(PVres_TotalAnnualCost(:));
    
else
    
    PVres_TotalAnnualCost = AnnuityFactor * (PV_SpecInvest + PVres_SpecInvest) * WPM2{m}(r,c)*1e-3 + (PVres_OnMCost * WPM2{m}(r,c)*1e-3);	% Total Annual Cost for 1 KW
    
end

%% PVres COE
if regcapcost==1
    for r=1:nr
        for c=1:nc
            if IRegion(r,c)==0 || IRegion(r,c)==27
                
                PVres_COE_Eff(r,c) = PVres_TotalAnnualCost_mean / PVres_TechPot_Cell_Eff(r,c); % $ / kWh
                PVres_COE_LF(r,c) = PVres_TotalAnnualCost_mean / PVres_TechPot_Cell_LF(r,c); % $ / kWh
                
                PVres_COE_Effu(r,c) = PVres_TotalAnnualCost_mean / PVres_TechPot_Cell_Eff(r,c); % $ / kWh
                PVres_COE_LFu(r,c) = PVres_TotalAnnualCost_mean / PVres_TechPot_Cell_LF(r,c); % $ / kWh
                
                PVres_COE_Effr(r,c) = PVres_TotalAnnualCost_mean / PVres_TechPot_Cell_Eff(r,c); % $ / kWh
                PVres_COE_LFr(r,c) = PVres_TotalAnnualCost_mean / PVres_TechPot_Cell_LF(r,c); % $ / kWh
                
            else
                
                PVres_COE_Eff(r,c) = PVres_TotalAnnualCost(IRegion(r,c)) / PVres_TechPot_Cell_Eff(r,c); % $ / kWh
                PVres_COE_LF(r,c) = PVres_TotalAnnualCost(IRegion(r,c)) / PVres_TechPot_Cell_LF(r,c); % $ / kWh
                
                PVres_COE_Effu(r,c) = PVres_TotalAnnualCost(IRegion(r,c)) / PVres_TechPot_Cell_Eff(r,c); % $ / kWh
                PVres_COE_LFu(r,c) = PVres_TotalAnnualCost(IRegion(r,c)) / PVres_TechPot_Cell_LF(r,c); % $ / kWh
                
                PVres_COE_Effr(r,c) = PVres_TotalAnnualCost(IRegion(r,c)) / PVres_TechPot_Cell_Eff(r,c); % $ / kWh
                PVres_COE_LFr(r,c) = PVres_TotalAnnualCost(IRegion(r,c)) / PVres_TechPot_Cell_LF(r,c); % $ / kWh
            end
        end
    end
    
else
    for r=1:nr
        for c=1:nc
            PVres_COE_Eff(r,c) = PVres_TotalAnnualCost / PVres_TechPot_Cell_Eff(r,c); % $ / kWh
            PVres_COE_LF(r,c) = PVres_TotalAnnualCost / PVres_TechPot_Cell_LF(r,c); % $ / kWh
            
            PVres_COE_Effu(r,c) = PVres_TotalAnnualCost / PVres_TechPot_Cell_Eff(r,c); % $ / kWh
            PVres_COE_LFu(r,c) = PVres_TotalAnnualCost / PVres_TechPot_Cell_LF(r,c); % $ / kWh
            
            PVres_COE_Effr(r,c) = PVres_TotalAnnualCost / PVres_TechPot_Cell_Eff(r,c); % $ / kWh
            PVres_COE_LFr(r,c) = PVres_TotalAnnualCost / PVres_TechPot_Cell_LF(r,c); % $ / kWh
        end
    end
end

PVres_Loadfac_Cell_prep = PV_Loadfac_Cell{13};
PVres_Loadfac_Cell_prepu = PV_Loadfac_Cell{13};
PVres_Loadfac_Cell_prepr = PV_Loadfac_Cell{13};

for r=1:nr
    for c=1:nc
        if AnnTechPotCellPVres(r,c)==0; %Total
            PVres_COE_LF(r,c) = NaN;
            PVres_COE_Eff(r,c) = NaN;
            PVres_Loadfac_Cell_prep(r,c) = NaN;
        end
        if AnnTechPotCellPVresu(r,c)==0 %Urban
            PVres_COE_LFu(r,c) = NaN;
            PVres_COE_Effu(r,c) = NaN;
            PVres_Loadfac_Cell_prepu(r,c) = NaN;
        end
        if AnnTechPotCellPVresr(r,c)==0 %Rural
            PVres_COE_LFr(r,c) = NaN;
            PVres_COE_Effr(r,c) = NaN;
            PVres_Loadfac_Cell_prepr(r,c) = NaN;
        end
    end
end

if setEconCeil==1
    PVres_COE_LF(PVres_COE_LF>CostCeil)=NaN;
    PVres_COE_Eff(PVres_COE_Eff>CostCeil)=NaN;
    
    PVres_COE_LFu(PVres_COE_LFu>CostCeil)=NaN;
    PVres_COE_Effu(PVres_COE_Effu>CostCeil)=NaN;
    
    PVres_COE_LFr(PVres_COE_LFr>CostCeil)=NaN;
    PVres_COE_Effr(PVres_COE_Effr>CostCeil)=NaN;
end

% figure(3);clf;imagesc(PVres_COE_Eff);axis image; colorbar
% figure(4);clf;imagesc(PVres_COE_LF);axis image; colorbar

%% PVres Cost curve prep IMAGE regions
clear SI
for i=1:26
        
    PVres_COE_IR{i} = PVres_COE_LF(IRind{i}); %Total
    PVres_Pot_IR{i} = AnnTechPotCellPVres(IRind{i});
    PVres_LF_IR{i} = PVres_Loadfac_Cell_prep(IRind{i});
    
    PVres_COE_IRu{i} = PVres_COE_LFu(IRind{i}); %Urban
    PVres_Pot_IRu{i} = AnnTechPotCellPVresu(IRind{i});
    PVres_LF_IRu{i} = PVres_Loadfac_Cell_prepu(IRind{i});
    
    PVres_COE_IRr{i} = PVres_COE_LFr(IRind{i}); %Rural
    PVres_Pot_IRr{i} = AnnTechPotCellPVresr(IRind{i});
    PVres_LF_IRr{i} = PVres_Loadfac_Cell_prepr(IRind{i});
end

for i=1:(nr*nc)
    PVres_COE_IR{27}(i) = PVres_COE_LF(i); %Total
    PVres_Pot_IR{27}(i) = AnnTechPotCellPVres(i);
    PVres_LF_IR{27}(i) = PVres_Loadfac_Cell_prep(i);
    
    PVres_COE_IRu{27}(i) = PVres_COE_LFu(i); %Urban
    PVres_Pot_IRu{27}(i) = AnnTechPotCellPVresu(i);
    PVres_LF_IRu{27}(i) = PVres_Loadfac_Cell_prepu(i);
    
    PVres_COE_IRr{27}(i) = PVres_COE_LFr(i); %Rural
    PVres_Pot_IRr{27}(i) = AnnTechPotCellPVresr(i);
    PVres_LF_IRr{27}(i) = PVres_Loadfac_Cell_prepr(i);
end

for i=1:26
    [PVres_COE_IR_S{i}, SI{i}] = sort(PVres_COE_IR{i}); %Total
    PVres_Pot_IR_S{i} = PVres_Pot_IR{i}(SI{i});
    PVres_LF_IR_S{i} = PVres_LF_IR{i}(SI{i});
    
    [PVres_COE_IR_Su{i}, SIu{i}] = sort(PVres_COE_IRu{i}); %Urban
    PVres_Pot_IR_Su{i} = PVres_Pot_IRu{i}(SIu{i});
    PVres_LF_IR_Su{i} = PVres_LF_IRu{i}(SIu{i});
    
    [PVres_COE_IR_Sr{i}, SIr{i}] = sort(PVres_COE_IRr{i}); %Rural
    PVres_Pot_IR_Sr{i} = PVres_Pot_IRr{i}(SIr{i});
    PVres_LF_IR_Sr{i} = PVres_LF_IRr{i}(SIr{i});
end

[PVres_COE_IR_S{27}, SI{27}] = sort(PVres_COE_IR{27}); %Total
PVres_Pot_IR_S{27} = PVres_Pot_IR{27}(SI{27});
PVres_LF_IR_S{27} = PVres_LF_IR{27}(SI{27});

[PVres_COE_IR_Su{27}, SIu{27}] = sort(PVres_COE_IRu{27}); %Urban
PVres_Pot_IR_Su{27} = PVres_Pot_IRu{27}(SIu{27});
PVres_LF_IR_Su{27} = PVres_LF_IRu{27}(SIu{27});

[PVres_COE_IR_Sr{27}, SIr{27}] = sort(PVres_COE_IRr{27}); %Rural
PVres_Pot_IR_Sr{27} = PVres_Pot_IRr{27}(SIr{27});
PVres_LF_IR_Sr{27} = PVres_LF_IRr{27}(SIr{27});

for i=1:27
    PVres_COE_IR_S{i} = PVres_COE_IR_S{i}(~isnan(PVres_COE_IR_S{i})); %Total
    PVres_Pot_IR_S{i} = PVres_Pot_IR_S{i}(~isnan(PVres_COE_IR_S{i}));
    PVres_LF_IR_S{i} = PVres_LF_IR_S{i}(~isnan(PVres_COE_IR_S{i}));
    
    PVres_COE_IR_Su{i} = PVres_COE_IR_Su{i}(~isnan(PVres_COE_IR_Su{i})); %Urban
    PVres_Pot_IR_Su{i} = PVres_Pot_IR_Su{i}(~isnan(PVres_COE_IR_Su{i}));
    PVres_LF_IR_Su{i} = PVres_LF_IR_Su{i}(~isnan(PVres_COE_IR_Su{i}));
    
    PVres_COE_IR_Sr{i} = PVres_COE_IR_Sr{i}(~isnan(PVres_COE_IR_Sr{i})); %Rural
    PVres_Pot_IR_Sr{i} = PVres_Pot_IR_Sr{i}(~isnan(PVres_COE_IR_Sr{i}));
    PVres_LF_IR_Sr{i} = PVres_LF_IR_Sr{i}(~isnan(PVres_COE_IR_Sr{i}));
    
end

for i=1:27
    RegEconPotPVres(i) = sum(PVres_Pot_IR_S{i}); % kWh / region / y
    RegEconPotPVresu(i) = sum(PVres_Pot_IR_Su{i}); % kWh / region / y
    RegEconPotPVresr(i) = sum(PVres_Pot_IR_Sr{i}); % kWh / region / y
end;


for i=1:27
    idx=find(~isnan(PVres_COE_IR{i}));
    PVres_COE_IR_wavg(i) = sum((PVres_COE_IR{i}(idx).*PVres_Pot_IR{i}(idx)))/sum(PVres_Pot_IR{i}(idx));
%    fprintf('%s %0.2f $/kWh\n',IMAGER{i},PVres_COE_IR_wavg(i))
end

%% PVres Cost curve prep Countries
clear SI
for i=1:(ncw-1)
        
    PVres_COE_C{i} = PVres_COE_LF(Cind{i}); %Total
    PVres_Pot_C{i} = AnnTechPotCellPVres(Cind{i});
    PVres_LF_C{i} = PVres_Loadfac_Cell_prep(Cind{i});
    
    PVres_COE_Cu{i} = PVres_COE_LFu(Cind{i}); %Urban
    PVres_Pot_Cu{i} = AnnTechPotCellPVresu(Cind{i});
    PVres_LF_Cu{i} = PVres_Loadfac_Cell_prepu(Cind{i});
    
    PVres_COE_Cr{i} = PVres_COE_LFr(Cind{i}); %Rural
    PVres_Pot_Cr{i} = AnnTechPotCellPVresr(Cind{i});
    PVres_LF_Cr{i} = PVres_Loadfac_Cell_prepr(Cind{i});
end

for i=1:(nr*nc)
    PVres_COE_C{ncw}(i) = PVres_COE_LF(i); %Total
    PVres_Pot_C{ncw}(i) = AnnTechPotCellPVres(i);
    PVres_LF_C{ncw}(i) = PVres_Loadfac_Cell_prep(i);
    
    PVres_COE_Cu{ncw}(i) = PVres_COE_LFu(i); %Urban
    PVres_Pot_Cu{ncw}(i) = AnnTechPotCellPVresu(i);
    PVres_LF_Cu{ncw}(i) = PVres_Loadfac_Cell_prepu(i);
    
    PVres_COE_Cr{ncw}(i) = PVres_COE_LFr(i); %Rural
    PVres_Pot_Cr{ncw}(i) = AnnTechPotCellPVresr(i);
    PVres_LF_Cr{ncw}(i) = PVres_Loadfac_Cell_prepr(i);
end

for i=1:(ncw-1)
    [PVres_COE_C_S{i}, SI{i}] = sort(PVres_COE_C{i}); %Total
    PVres_Pot_C_S{i} = PVres_Pot_C{i}(SI{i});
    PVres_LF_C_S{i} = PVres_LF_C{i}(SI{i});
    
    [PVres_COE_C_Su{i}, SIu{i}] = sort(PVres_COE_Cu{i}); %Urban
    PVres_Pot_C_Su{i} = PVres_Pot_Cu{i}(SIu{i});
    PVres_LF_C_Su{i} = PVres_LF_Cu{i}(SIu{i});
    
    [PVres_COE_C_Sr{i}, SIr{i}] = sort(PVres_COE_Cr{i}); %Rural
    PVres_Pot_C_Sr{i} = PVres_Pot_Cr{i}(SIr{i});
    PVres_LF_C_Sr{i} = PVres_LF_Cr{i}(SIr{i});
end

[PVres_COE_C_S{ncw}, SI{ncw}] = sort(PVres_COE_C{ncw}); %Total
PVres_Pot_C_S{ncw} = PVres_Pot_C{ncw}(SI{ncw});
PVres_LF_C_S{ncw} = PVres_LF_C{ncw}(SI{ncw});

[PVres_COE_C_Su{ncw}, SIu{ncw}] = sort(PVres_COE_Cu{ncw}); %Urban
PVres_Pot_C_Su{ncw} = PVres_Pot_Cu{ncw}(SIu{ncw});
PVres_LF_C_Su{ncw} = PVres_LF_Cu{ncw}(SIu{ncw});

[PVres_COE_C_Sr{ncw}, SIr{ncw}] = sort(PVres_COE_Cr{ncw}); %Rural
PVres_Pot_C_Sr{ncw} = PVres_Pot_Cr{ncw}(SIr{ncw});
PVres_LF_C_Sr{ncw} = PVres_LF_Cr{ncw}(SIr{ncw});

for i=1:ncw
    PVres_COE_C_S{i} = PVres_COE_C_S{i}(~isnan(PVres_COE_C_S{i})); %Total
    PVres_Pot_C_S{i} = PVres_Pot_C_S{i}(~isnan(PVres_COE_C_S{i}));
    PVres_LF_C_S{i} = PVres_LF_C_S{i}(~isnan(PVres_COE_C_S{i}));
    
    PVres_COE_C_Su{i} = PVres_COE_C_Su{i}(~isnan(PVres_COE_C_Su{i})); %Urban
    PVres_Pot_C_Su{i} = PVres_Pot_C_Su{i}(~isnan(PVres_COE_C_Su{i}));
    PVres_LF_C_Su{i} = PVres_LF_C_Su{i}(~isnan(PVres_COE_C_Su{i}));
    
    PVres_COE_C_Sr{i} = PVres_COE_C_Sr{i}(~isnan(PVres_COE_C_Sr{i})); %Rural
    PVres_Pot_C_Sr{i} = PVres_Pot_C_Sr{i}(~isnan(PVres_COE_C_Sr{i}));
    PVres_LF_C_Sr{i} = PVres_LF_C_Sr{i}(~isnan(PVres_COE_C_Sr{i}));
    
end

for i=1:ncw
    CEconPotPVres(i) = sum(PVres_Pot_C_S{i}); % kWh / country / y
    CEconPotPVresu(i) = sum(PVres_Pot_C_Su{i}); % kWh / country/ y
    CEconPotPVresr(i) = sum(PVres_Pot_C_Sr{i}); % kWh / country / y
end;

%% PVres Cost curve prep for paper graph
clear SI

%Africa
PaperReg{1} = [8,9,10,16];
%Europe
PaperReg{2} = [11,12,13];
%India
PaperReg{3} = [18,25];
%SAM
PaperReg{4} = [3,4,5,6];
%MENA
PaperReg{5} = [7,17];
%NAM
PaperReg{6} = [1,2];
%POECD
PaperReg{7} = [23,24];
%Russia
PaperReg{8} = [16,14,15];
%RAsia
PaperReg{9} = [22,19,21];
%China
PaperReg{10} = 20;

%Find cell per Paper reg
for i=1:10
    
    for j=1:numel(PaperReg{i})
        PRind2{i}{j} = find(IRegion(:)==PaperReg{i}(j)); %PVres Paper regions
    end
    PRind{i} = vertcat(PRind2{i}{:});
    
end

%% sum to 10 regions
clear SI
for i=1:10
    PVres_COE_IR_p{i} = PVres_COE_LF(PRind{i}); %Total
    PVres_Pot_IR_p{i} = AnnTechPotCellPVres(PRind{i});
    PVres_LF_IR_p{i} = PVres_Loadfac_Cell_prep(PRind{i});
end

for i=1:(nr*nc)
    PVres_COE_IR_p{11}(i) = PVres_COE_LF(i); %Total
    PVres_Pot_IR_p{11}(i) = AnnTechPotCellPVres(i);
    PVres_LF_IR_p{11}(i) = PVres_Loadfac_Cell_prep(i);
end

for i=1:10
    [PVres_COE_IR_S_p{i}, SI{i}] = sort(PVres_COE_IR_p{i}); %Total
    PVres_Pot_IR_S_p{i} = PVres_Pot_IR_p{i}(SI{i});
    PVres_LF_IR_S_p{i} = PVres_LF_IR_p{i}(SI{i});
end

[PVres_COE_IR_S_p{11}, SI{11}] = sort(PVres_COE_IR_p{11}); %Total
PVres_Pot_IR_S_p{11} = PVres_Pot_IR_p{11}(SI{11});
PVres_LF_IR_S_p{11} = PVres_LF_IR_p{11}(SI{11});

for i=1:11
    PVres_COE_IR_S_p{i} = PVres_COE_IR_S_p{i}(~isnan(PVres_COE_IR_S_p{i})); %Total
    PVres_Pot_IR_S_p{i} = PVres_Pot_IR_S_p{i}(~isnan(PVres_Pot_IR_S_p{i}));
    PVres_LF_IR_S_p{i} = PVres_LF_IR_S_p{i}(~isnan(PVres_LF_IR_S_p{i}));
end

%% Weighted average COE per region
i=1;
for i=1:10
    idx=find(~isnan(PVres_COE_IR_p{i}));
    PVres_COE_IR_p_wavg(i) = sum((PVres_COE_IR_p{i}(idx).*PVres_Pot_IR_p{i}(idx)))/sum(PVres_Pot_IR_p{i}(idx));
    %fprintf('%0.2f $/kWh\n',PVres_COE_IR_p_wavg(i))
end

fname = sprintf('%s\\COEwa.xlsx',root);
%xlswrite(fname,PVres_COE_IR_p_wavg')
