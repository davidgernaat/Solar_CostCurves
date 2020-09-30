%% Suitability factors

% SuitFact1 = 0.01,
% SuitFact2 = 0.05,
% SuitFact3 = 0.10,
% SuitFact4 = 0.20,
% SuitFact5 = 0.25,
% SuitFact6 = 0.30,
% SuitFact10 = 1.00;

% !1	agricultural land Hoogwijk: 0.01
% !2	extensive grassland 0.05
% !3	carbon plantations
% !4	regrowth forest (Abandoning)
% !5	regrowth forest (Timber)
% !6	biofuel
% !7	ice
% !8	tundra 0.01
% !9	wooded tundra
% !10	boreal forest
% !11	cool conifer
% !12	temp. mixed forest
% !13	temp.deciduous forest
% !14	warm mixed
% !15	grassland/steppe 0.01
% !16	hot desert 0.05
% !17	scrubland 0.01
% !18	savanna 0.01
% !19	tropical woodland
% !20	tropical forest

for r=1:nr
    for c=1:nc
        if GLCT(r,c)==7 || GLCT(r,c)==9 || GLCT(r,c)==10 || GLCT(r,c)==11 || GLCT(r,c)==12 || GLCT(r,c)==13 || GLCT(r,c)==14 || GLCT(r,c)==20
            SuitabilityFactor(r,c) = 0;
        elseif GLCT(r,c)==1
            SuitabilityFactor(r,c) = 0.01;
        elseif GLCT(r,c)==2
            SuitabilityFactor(r,c) = 0.05;
        elseif GLCT(r,c)==4
            SuitabilityFactor(r,c) = 0;
        elseif GLCT(r,c)==5
            SuitabilityFactor(r,c) = 0;
        elseif GLCT(r,c)==6
            SuitabilityFactor(r,c) = 0;
        elseif GLCT(r,c)==8
            SuitabilityFactor(r,c) = 0.01;
        elseif GLCT(r,c)==15
            SuitabilityFactor(r,c) = 0.01;
        elseif GLCT(r,c)==16
            SuitabilityFactor(r,c) = 0.05;
        elseif GLCT(r,c)==17
            SuitabilityFactor(r,c) = 0.01;
        elseif GLCT(r,c)==18
            SuitabilityFactor(r,c) = 0.01;
        elseif GLCT(r,c)==19
            SuitabilityFactor(r,c) = 0;
        else
            SuitabilityFactor(r,c) = 0;
        end
    end
end

%% BIOres if zero than oke
for r=1:nr
    for c=1:nc
        if BIOres(r,c)==0
            BioMult(r,c)=1;
        else
            BioMult(r,c)=0;
        end
    end
end

%% Builtup
for r=1:nr
    for c=1:nc
        BuildupInv(r,c)=min(1, max(0, 1 - Buildup(r,c)));

        if BuildupInv(r,c) < 0.9
            BuildupInv(r,c) = 0;
        end
    end
end

% figure(1);clf;imagesc(BuildupInv); colorbar

% Fraction of cell area with slope <3%

%% ExclFactor
for r=1:nr
    for c=1:nc
        ExclFactor(r,c) = SuitabilityFactor(r,c) * BioMult(r,c) * BuildupInv(r,c) * Slope(r,c);
    end
end
% figure(2);clf;imagesc(ExclFactor); colormap(jet); colorbar

% Exclusion factor map
% txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.2f\nNODATA_value\t%d',720,360,-180,-90,0.5,-99);
% file = fullfile(root, sprintf('\\output\\Suitability_map.asc'));
% dlmwrite(file,txt,'');
% dlmwrite(file,ExclFactor,'-append','delimiter',' '); % 0-1 suitability factor

%% CSP Geographic potential
m2tokm2conv = 1000000;		%1 km2 = 1.000.000 m2
for m=13
    for r=1:nr
        for c=1:nc
            CSP_GeoPot_Cell{m}(r,c) = Avg_DNI_NASA_adj{m}(r,c) * Area(r,c) * m2tokm2conv * 365 ...
                * ExclFactor(r,c); % kWh / cell / y
        end
    end
end

for i=1:26
    RegGeoPotCSP(i) = sum(CSP_GeoPot_Cell{13}(IRind{i})); % kWh / region / y
end
RegGeoPotCSP(27) =sum(RegGeoPotCSP(1:26));

for i=1:numel(ISOGDP(:,1))
    CGeoPotCSP(i) = sum(CSP_GeoPot_Cell{13}(Cind{i})); % kWh / country / y
end
dt = sum(CGeoPotCSP(1:end));
CGeoPotCSP(end+1) = dt;

%% PV Geographical potential
for m=13
    for r=1:nr
        for c=1:nc
            PV_GeoPot_Cell{m}(r,c) = Gl_Horiz_NASA{m}(r,c) * Area(r,c) * m2tokm2conv * 365 ...
                * ExclFactor(r,c); %kWh / cell / y
        end
    end
end

for i=1:26
    RegGeoPotPV(i) = sum(PV_GeoPot_Cell{m}(IRind{i})); % kWh / region / y
end
RegGeoPotPV(27) =sum(RegGeoPotPV(1:26));

for i=1:numel(ISOGDP(:,1))
    CGeoPotPV(i) = sum(PV_GeoPot_Cell{m}(Cind{i})); % kWh / country / y
end
dt = sum(CGeoPotPV(1:end));
CGeoPotPV(end+1) = dt;

%% PV res
%Formula from Tiene's work

% if regbeta==0 %Switched off now, urban/rural numbers are not used anymore, global ones instead
%     
%     for i=1:26
%         Bu(i)=0.5;
%         Br(i)=7;
%     end
%     
% else
%     Bu = [0.4 0.4 0.5 0.5 0.5 0.5 0.57 0.6 0.6 0.5 0.46 0.73 0.2 0.73 0.6 0.6 0.57 0.96 0.2 0.5 0.96 0.96 0.2 0.8 0.96 0.5 0.5];
%     Br = [0.7 0.7 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.7 0.72 0.9 0.5 0.9 1 1 0.8 1 0.5 0.7 1 1 0.5 0.9 1 0.7 0.7];
%     
% end

for R=1:27
    %FloorArea_Region(R)	=  resFloorspace_TURQ(45,R,2) * Pop_TURQ(45,R,2) + resFloorspace_TURQ(45,R,3) * Pop_TURQ(45,R,3);%Urban/rural
    FloorArea_Region(R)	=  resFloorspace_TURQ(45,R,1) * Pop_TURQ(45,R,1); %total
end

%fprintf('Global floor space: %.0f billion km2\n',sum(FloorArea_Region(:))*1e-9)

for R=1:27
    %RoofArea_Region(R)	=  Bu(R) * resFloorspace_TURQ(45,R,2) * Pop_TURQ(45,R,2) + Br(R) * resFloorspace_TURQ(45,R,3) * Pop_TURQ(45,R,3);
    RoofArea_Region(R)	=  Beta_fs_ra(R) * resFloorspace_TURQ(45,R,1) * Pop_TURQ(45,R,1); % total
end

% fprintf('Total Roof area global: %.0f km2\n',(sum(RoofArea_Region(1:26))*1e-6))
% fprintf('suitable Roof area global: %.0f km2\n',(sum(RoofArea_Region(1:26))*1e-6*0.32))
% fprintf('USA: %.0f km2\n',(RoofArea_Region(2))*1e-6*0.32)
% fprintf('EU: %.0f km2\n',(RoofArea_Region(11)+RoofArea_Region(12))*1e-6*0.32)
% fprintf('Brazil : %.0f km2\n',(RoofArea_Region(5))*1e-6*0.32)

%% Put Floor area in map
for j=1:3
    FAmap{j} = zeros(nr,nc);
    for r=1:nr
        for c=1:nc
            if IRegion(r,c)==0; continue; end
            if j==1
                FAmap{j}(r,c) = resFloorspace_TURQ(45,IRegion(r,c),1); % m2 / capita total
            elseif j==2
                FAmap{j}(r,c) = resFloorspace_TURQ(45,IRegion(r,c),2); % m2 / capita urban
            elseif j==3
                FAmap{j}(r,c) = resFloorspace_TURQ(45,IRegion(r,c),3); % m2 / capita rural
            end
        end
    end
end
% imagesc(FAmap{1}); axis image

%% Scale urbpop and rurpop so that it matches the total of Pop_Turq
for R=1:26;
    totpop(IRind{R}) = totpop(IRind{R})*(Pop_TURQ(45,R,1)/sum(totpop(IRind{R})));
    urbpop(IRind{R}) = urbpop(IRind{R})*(Pop_TURQ(45,R,2)/sum(urbpop(IRind{R})));
    rurpop(IRind{R}) = rurpop(IRind{R})*(Pop_TURQ(45,R,3)/sum(rurpop(IRind{R})));
end

%%
for i=1:26
    TotPop(i) = sum(totpop(find(IRegion(:)==i))); % Total number of people per region
    TotUrb(i) = sum(urbpop(find(IRegion(:)==i))); % Urban number of people per region
    TotRur(i) = sum(rurpop(find(IRegion(:)==i))); % Rural number of people per region
end
TotPop(27) =sum(TotPop(1:26));
TotUrb(27) =sum(TotUrb(1:26));
TotRur(27) =sum(TotRur(1:26));

for i=1:numel(ISOGDP(:,1))
    CTotPop(i) = sum(totpop(find(IRegion(:)==ISOGDP(i,1)))); % Total number of people per country
    CTotUrb(i) = sum(urbpop(find(IRegion(:)==ISOGDP(i,1)))); % Urban number of people per country
    CTotRur(i) = sum(rurpop(find(IRegion(:)==ISOGDP(i,1)))); % Rural number of people per country
end
dt = sum(CTotPop(1:end));
CTotPop(end+1) = dt;
dt = sum(CTotUrb(1:end));
CTotUrb(end+1) = dt;
dt = sum(CTotRur(1:end));
CTotRur(end+1) = dt;

%% Same formula from Tiene
RAmap = zeros(nr,nc);
RAmapu = zeros(nr,nc);
RAmapr = zeros(nr,nc);

for r=1:nr
    for c=1:nc
        if (IRegion(r,c))==0
            RAmap(r,c) = 0;
            RAmapu(r,c) = 0;
            RAmapr(r,c) = 0;
        else
            RAmap(r,c)  = Beta_fs_ra(IRegion(r,c)) * totpop(r,c) * FAmap{1}(r,c); % Total m2 roof area per cell
            RAmapu(r,c) = Beta_fs_ra(IRegion(r,c)) * urbpop(r,c) * FAmap{2}(r,c); % Total m2 roof area per cell urban
            RAmapr(r,c) = Beta_fs_ra(IRegion(r,c)) * rurpop(r,c) * FAmap{3}(r,c); % Total m2 roof area per cell rural
        end
    end
end

fprintf('Global %.0f km2\n',(sum(RAmap(:)))*1e-6*0.32);
fprintf('Roof  US %.0f km2\n',(sum(RAmap(IRind{2})))*1e-6*0.32);
fprintf('EU %.0f km2\n',(sum(RAmap(IRind{11}))+sum(RAmap(IRind{12})))*1e-6*0.32);
fprintf('WEU %.0f km2\n',(sum(RAmap(IRind{11})))*1e-6*0.32);
fprintf('Brazil %.0f km2\n',(sum(RAmap(IRind{5})))*1e-6*0.32);
fprintf('Brazil %.0f km2\n',(sum(RAmap(IRind{5})))*1e-6);
fprintf('Switserland %.0f km2\n',(sum(RAmap(Cind{171})))*1e-6);
fprintf('Spain %.0f km2\n',(sum(RAmap(Cind{166})))*1e-6);

% Roof area map
% txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.2f\nNODATA_value\t%d',720,360,-180,-90,0.5,-99);
% file = fullfile(root, sprintf('\\output\\Rooftop_map (m2 per cell).asc'));
% dlmwrite(file,txt,'');
% dlmwrite(file,RAmap,'-append','delimiter',' '); % 0-1 suitability factor

%% PVres TheoPot
for m=13
    for r=1:nr
        for c=1:nc
            PVres_TheoPot{m}(r,c)  = RAmap(r,c)  * Gl_Horiz_NASA{m}(r,c) * 365;   % kWh / y
            PVres_TheoPotu{m}(r,c) = RAmapu(r,c) * Gl_Horiz_NASA{m}(r,c) * 365; % kWh / y urban
            PVres_TheoPotr{m}(r,c) = RAmapr(r,c) * Gl_Horiz_NASA{m}(r,c) * 365; % kWh / y rural
        end
    end
end
for i=1:26
    RegTheoPotPVres(i)  = sum(PVres_TheoPot{13}(IRind{i}));   % kWh / region / y
    RegTheoPotPVresu(i) = sum(PVres_TheoPotu{13}(IRind{i})); % kWh / region / y urban
    RegTheoPotPVresr(i) = sum(PVres_TheoPotr{13}(IRind{i})); % kWh / region / y rural
end
RegTheoPotPVres(27)  =sum(RegTheoPotPVres(1:26));
RegTheoPotPVresu(27) =sum(RegTheoPotPVres(1:26));
RegTheoPotPVresr(27) =sum(RegTheoPotPVres(1:26));

for i=1:numel(ISOGDP(:,1))
    CTheoPotPVres(i)  = sum(PVres_TheoPot{13}(Cind{i}));   % kWh / country / y
    CTheoPotPVresu(i) = sum(PVres_TheoPotu{13}(Cind{i})); % kWh / country / y urban
    CTheoPotPVresr(i) = sum(PVres_TheoPotr{13}(Cind{i})); % kWh / country / y rural
end
CTheoPotPVres(end+1) =sum(CTheoPotPVres(1:end));
CTheoPotPVresu(end+1) =sum(CTheoPotPVresu(1:end));
CTheoPotPVresr(end+1) =sum(CTheoPotPVresr(1:end));

%% PVres GeoPot
PVres_RoofUseFactor = 0.32; %Architectual obstrusian + bad roofs

for m=13
    for r=1:nr
        for c=1:nc
            PVres_GeoPot{m}(r,c)  = RAmap(r,c)  * Gl_Horiz_NASA{m}(r,c) * 365 * PVres_RoofUseFactor;   % kWh / y
            PVres_GeoPotu{m}(r,c) = RAmapu(r,c) * Gl_Horiz_NASA{m}(r,c) * 365 * PVres_RoofUseFactor; % kWh / y urban
            PVres_GeoPotr{m}(r,c) = RAmapr(r,c) * Gl_Horiz_NASA{m}(r,c) * 365 * PVres_RoofUseFactor; % kWh / y rural
        end
    end
end

for i=1:26
    RegGeoPotPVres(i)  = sum(PVres_GeoPot{13}(IRind{i}));   % kWh / region / y
    RegGeoPotPVresu(i) = sum(PVres_GeoPotu{13}(IRind{i})); % kWh / region / y urban
    RegGeoPotPVresr(i) = sum(PVres_GeoPotr{13}(IRind{i})); % kWh / region / y rural
end
RegGeoPotPVres(27)  = sum(RegGeoPotPVres(1:26));
RegGeoPotPVresu(27) = sum(RegGeoPotPVres(1:26));
RegGeoPotPVresr(27) = sum(RegGeoPotPVres(1:26));

for i=1:numel(ISOGDP(:,1))
    CGeoPotPVres(i)  = sum(PVres_GeoPot{13}(Cind{i}));   % kWh / country / y
    CGeoPotPVresu(i) = sum(PVres_GeoPotu{13}(Cind{i})); % kWh / country / y urban
    CGeoPotPVresr(i) = sum(PVres_GeoPotr{13}(Cind{i})); % kWh / country / y rural
end
CGeoPotPVres(end+1)  = sum(CGeoPotPVres(1:end));
CGeoPotPVresu(end+1) = sum(CGeoPotPVresu(1:end));
CGeoPotPVresr(end+1) = sum(CGeoPotPVresr(1:end));

%%
% win=0;
% a = EEZ==11;
% cols = find(sum(a)>0);
% rows = find(sum(a')>0);
% fc = cols(1) - win;
% lc = cols(end) + win;
% fr = rows(1) - win;
% lr = rows(end) + win;
% 
% b=TechPotCell{13};
% b(find(~a))=-1;
% 
% bw = b(fr:lr,fc:lc);
% 
% figure(1);clf;imagesc(bw);colorbar; axis image
