%% CSP PV and PVresidential cost supply curve generator
clear all

root = 'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER\CSP_PV_PVres';

fname = sprintf('%s\\data\\input_data.mat', root);
load(fname)

[nr nc] = size(IRegion);

%% Settings
% 0=no/1=yes
setoutput=1;                % Do you want output (cost curves, max prod, load curves)?
outputmaps=0;               % Do you want to output ASCII maps?
outputcountrypots=0;        % Do you want to output country potentials?
setCSPPL=1;                 % Do you want to include CSP powerline costs?
setPVPL=1;                  % Do you want to include PV powerline costs?
setEconCeil=0;              % Do you want to add a cost ceiling for economic potential?
multifactor=0;              % Do you want to use multiple powerline setting?
regcapcost=1;               % Do you want to include regional capital cost (1) or a global uniform capital cost (0) 
regbeta=1;                  % Do you want to include regional beta to translate floorspace in roof area (1) or global betas (0)
CostCeil=0.15;               % Level of cost ceiling $ / kWh

%% Literature numbers
fprintf('\nGlobal    Theoretical Potential PV Hoogwijk %.2f PWh', 176e3);
fprintf('\nGlobal     Geographic Potential PV Hoogwijk %.2f PWh', 3.5e3);
fprintf('\nGlobal  Geographic Potential PVres Hoogwijk %.2f PWh', 0.6e2);
fprintf('\nGlobal      Technical Potential PV Hoogwijk %.2f PWh', 366);
fprintf('\nGlobal   Technical Potential PVres Hoogwijk %.2f PWh\n\n', 6);

fprintf('Global        Technical Potential PV REMIND %.2f PWh', 3370 * 0.277777778);
fprintf('\nGlobal       Technical Potential CSP REMIND %.2f PWh\n\n', 2362 * 0.277777778);

%% Theoretical potential
disp('Theoretical potential')
Solar_TheoPot
% Method:      	!DNI cutoff for CSP is 3 kwh/m2/day; PV works with wider range and can also use diffuse radiation
% Read in DNI from NASA dataset, then adjust for cloudy days (NASA NO_SUN dataset)
% For each cell check if DNI is below CSP cutoff, if yes than set to PV only
% If Adj_DNI<3 then cell good only for PV (TechMarker=1)
% If 3<Adj_DNI<5 then CSP-PV competition exists
% If Adj_DNI>5 then preferential CSP use

fprintf('\nGlobal Theoretical Potential CSP %.2f PWh\n', sum(AnnPotCellDNI(:))*1e-12);
fprintf('Global Theoretical Potential PV  %.2f PWh\n\n', sum(AnnPotCellPV(:))*1e-12);

%% Geographical potential
disp('Geographical potential')
Solar_GeoPot
% Method:      	Theoretical potential adjusted for suitability of land use type.
% Land use class suitability excludes urban areas, bioreserves, land with more than 3% slope
% and different fractions of Land Cover Type is allowed (deserts, grasslands, savannas, etc)

fprintf('\nGlobal Geographic Potential CSP    %.2f PWh\n', sum(CSP_GeoPot_Cell{13}(:))*1e-12);
fprintf('Global Geographic Potential PV     %.2f PWh\n', sum(PV_GeoPot_Cell{13}(:))*1e-12);
fprintf('Global Geographic Potential POP PVres  %.2f PWh\n', sum(PVres_GeoPot{13}(:))*1e-12);
fprintf('Global Geographic Potential POP PVres Urban %.2f PWh\n', sum(PVres_GeoPotu{13}(:))*1e-12);
fprintf('Global Geographic Potential POP PVres Rural %.2f PWh\n', sum(PVres_GeoPotr{13}(:))*1e-12);
fprintf('Roof area total %.2f km2\n', sum(RAmap(:))*1e-6);
fprintf('Roof area suitable %.2f km2\n', sum(RAmap(:))*1e-6*PVres_RoofUseFactor);
fprintf('Roof area urban %.2f km2\n\n', sum(RAmapu(:))*1e-6*PVres_RoofUseFactor);

%% Technical potential
disp('Technical potential')
Solar_TechPot
% Method:      	geographhical potential multiplied by land use factor and the solar-to-electric efficiency of each technology.
% references:   Koberle 2013 Thesis
% Trieb et al 2009
% Peters et al 2011: Shedding light on solar technologies—A techno-economic assessment and its policy implications (in Energy Policy)
% Marion et al 2005: Performance Parameters for Grid-Connected PV Systems (NREL - 31st IEEE Photovoltaics Specialists, Conference and Exhibition Lake Buena Vista, Florida, January 3-7, 2005)
% Ong et al 2013: Land-Use Requirements for Solar Power Plants in the United States (NREL Technical Report NREL/TP-6A20-56290 June 2013

fprintf('\nGlobal Technical Potential CSP   %.2f PWh\n', sum(AnnTechPotCellCSP(:))*1e-12);
fprintf('Global Technical Potential PV    %.2f PWh\n', sum(AnnTechPotCellPV(:))*1e-12);
fprintf('Global Technical Potential PVres %.2f PWh\n', sum(AnnTechPotCellPVres(:))*1e-12);
fprintf('Global Technical Potential PVres Urban %.2f PWh\n', sum(AnnTechPotCellPVresu(:))*1e-12);
fprintf('Global Technical Potential PVres Rural %.2f PWh\n\n', sum(AnnTechPotCellPVresr(:))*1e-12);

%% Economic potential
disp('Economic potential')
Solar_EconPot

fprintf('\nGlobal Economic <%.2f $/kWh Potential CSP   %.2f PWh\n', CostCeil, sum(CSP_Pot_IR_S{end}(:))*1e-12);
fprintf('Global Economic <%.2f $/kWh Potential PV    %.2f PWh\n', CostCeil,sum(PV_Pot_IR_S{end}(:))*1e-12);
fprintf('Global Economic <%.2f $/kWh Potential PVres %.2f PWh\n', CostCeil,sum(PVres_Pot_IR_S{end}(:))*1e-12);
fprintf('Global Economic <%.2f $/kWh Potential PVres Urban %.2f PWh\n', CostCeil,sum(PVres_Pot_IR_Su{end}(:))*1e-12);
fprintf('Global Economic <%.2f $/kWh Potential PVres Rural %.2f PWh\n\n', CostCeil,sum(PVres_Pot_IR_Sr{end}(:))*1e-12);

%% Overview
disp('Show')
% Solar_Overview_CSP
% Solar_Overview_PV
% maps

%% Costsupply
disp('Cost Supply')
Solar_CostSupply

%% Output maps
if outputmaps
    disp('Generating output maps')
    Solar_output_maps
end

%% Output country potentials
if outputcountrypots
    disp('Country potentials output')
    Solar_PotsPerCountry
end