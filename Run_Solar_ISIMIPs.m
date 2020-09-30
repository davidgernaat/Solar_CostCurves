%% Run ISIMIPs
clear all

root = 'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER\CSP_PV_PVres';
root2 = 'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER';

root_rsds = sprintf('%s\\data\\ISIMIP\\Irradiance',root);
root_tas = sprintf('%s\\data\\ISIMIP\\Temperature',root);
root_ws = sprintf('%s\\data\\ISIMIP\\Wind',root);

% 1. rsds files : are in W m-2 (irradiance)
% 2. tas files : are in Kelvin, so 273.15 must be subtracted from the values to get degree Celsius (K-  273.15=oC)
% 3. ws files: are in m s-1 (wind speeds)

fname=sprintf('%s\\data\\imagemask_land.mat',root);
load(fname)

%% File list
disp('Read file list')
fnames_irr = dir(root_rsds);
fnames_tas = dir(root_tas);
fnames_ws = dir(root_ws);

c=0;
for j=1:numel(fnames_irr)
    if j<3; continue; end;
    c=c+1;
    pathnames_irr{c} = sprintf('%s\\%s',root_rsds,(fnames_irr(j).name));
    pathnames_tas{c} = sprintf('%s\\%s',root_tas,(fnames_tas(j).name));
    pathnames_ws{c} = sprintf('%s\\%s',root_ws,(fnames_ws(j).name));
end

c=0;
% clear names
% clc
for j=3:numel(fnames_irr)
    c=c+1;
    
    if c<6
        %         fprintf('%s\n',fnames_irr(j).name(10:64))
        names{c} = fnames_irr(j).name(10:64);
    elseif c>5 && c<12
        %         fprintf('%s\n',fnames_irr(j).name(10:64))
        names{c} = fnames_irr(j).name(10:64);
    elseif c>11 && c<16
        %         fprintf('%s\n',fnames_irr(j).name(10:64))
        names{c} = fnames_irr(j).name(10:64);
    elseif c>15
        %         fprintf('%s\n',fnames_irr(j).name(10:64))
        names{c} = fnames_irr(j).name(10:64);
    end
    
    names_split{c} =strsplit(names{c},'_');
    
    %fprintf('%s %s %s\n',names_split{c}{1},names_split{c}{2},names_split{c}{5})
end

for i=1:numel(names_split)
    GCMID{i}=names_split{i}{1};
    RCPID{i}=names_split{i}{2};
    TIMEID{i}=names_split{i}{5};
end

%% step 1: read netcdf
% i=6; %hadgem2-es hist
% i=10; %hadgem2-es rcp60 2100

for i=1:numel(names)
    fprintf('Reading file #%d of %d\n',i, numel(names))
    clear CM CMc dt
    
    %% Read, reformat and convert irradiance
    %% SOLAR
    disp('Solar data prep')
    %finfo = ncinfo(pathnames{i});
    %finfo.Variables.Name
    
    CM_data = ncread(pathnames_irr{i},'rsds');
    
    % step 2: convert to correct format
    for m=1:12
        for r=1:360
            for c=1:720
                CM{m}(r,c) = CM_data(c,r,m);
            end
        end
    end
    
    for r=1:360
        for c=1:720
            for m=1:12
                dt(m) = CM{m}(r,c);
            end
            CM{13}(r,c) = mean(dt);
        end
    end
    
    % step 3: converting units
    % kwh/m2/day --> W/m2 conversion
    % 1 kWh/m^2/day
    % =
    % 1 kW*hour/m^2/day
    % =
    % 1000 W*hour/m^2/day
    % =
    % 1000 W/m^2 * hour/day
    % =
    % 1000 W/m^2 * hour/(24 hours)
    % =
    % 1000 W/m^2 * (1 / (24))
    % =
    % 1000 W/m^2 * (.041667)
    % =
    % 41.667 W/m^2
    
    for m=1:13
        for r=1:360
            for c=1:720
                CMc{m}(r,c) = CM{m}(r,c) / 41.667;
            end
        end
    end
    
    % figure(1);clf;imagesc(CMc{13});colorbar
    
    %% TEMPERATURE
    disp('Temperature data prep')
    %finfo = ncinfo(pathnames_tas{i});
    %finfo.Variables.Name
    
    tas_data = ncread(pathnames_tas{i},'tas');
    
    % step 2: convert to correct format
    for m=1:12
        for r=1:360
            for c=1:720
                tas{m}(r,c) = tas_data(c,r,m);
            end
        end
    end
    
    for r=1:360
        for c=1:720
            for m=1:12
                dt(m) = tas{m}(r,c);
            end
            tas{13}(r,c) = mean(dt);
        end
    end
    
    % step 3: converting units
    for m=1:13
        for r=1:360
            for c=1:720
                tasc{m}(r,c) = tas{m}(r,c) - 273.15; %Kelvin --> Celcius
            end
        end
    end
    
    %figure(1);clf;imagesc(tasc{13});colorbar
    
    %% Wind Speeds
    disp('Wind data prep')
    %finfo = ncinfo(pathnames_ws{i});
    %finfo.Variables.Name
    
    ws_data = ncread(pathnames_ws{i},'sfcWind');
    
    % step 2: convert to correct format
    for m=1:12
        for r=1:360
            for c=1:720
                ws{m}(r,c) = ws_data(c,r,m);
            end
        end
    end
    
    for r=1:360
        for c=1:720
            for m=1:12
                dt(m) = ws{m}(r,c);
            end
            ws{13}(r,c) = mean(dt);
        end
    end
    
    % step 3: converting units
    for m=1:13
        for r=1:360
            for c=1:720
                wsc{m}(r,c) = ws{m}(r,c); %m/s
            end
        end
    end
    
    %figure(1);clf;imagesc(wsc{13});colorbar
    
    %% step 4: running Solar_cc with Climate model data
    disp('Running Cost-supply')
    
    [CostCurveSmthPV{i}, CostCurveSmthCSP{i}, CostCurveSmthPVres{i}, ...
        MaxProdPV{i}, MaxProdCSP{i}, MaxProdPVres{i}, ...
        LFCurveSmthPV{i}, LFCurveSmthCSP{i}, LFCurveSmthPVres{i}, ...
        AnnPotCellDNI{i}, CSP_GeoPot_Cell{i}, AnnTechPotCellCSP{i}, CSP_COE_PL{i}, CSP_COE_LF{i}, ...
        AnnPotCellPV{i}, PV_GeoPot_Cell{i}, AnnTechPotCellPV{i}, PV_COE_PLh{i}, PV_COE_LF{i}, ...
        PVres_TheoPot{i}, PVres_GeoPot{i}, AnnTechPotCellPVres{i}, PVres_COE_LF{i} ] = Solar_cc_ISIMIP(root,CMc,CMc,tasc,wsc);
    
%     a=PV_COE_LF{i};
%     a(a>0.3)=NaN;
%     figure(3);clf;imagesc(a);axis image;colorbar
    
    %% step 5: Deal with IAM output
    
    % POLES
    C2R_fname = sprintf('%s\\data\\ISIMIP\\Modelregionallocation\\Country_to_POLESregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\POLES\\Solar\\Solar_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellCSP{i},CSP_COE_PL{i},GCMID{i},RCPID{i},TIMEID{i},'CSP')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPV{i},PV_COE_PLh{i},GCMID{i},RCPID{i},TIMEID{i},'PV')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPVres{i},PVres_COE_LF{i},GCMID{i},RCPID{i},TIMEID{i},'PVres')
    
    % REMIND
    C2R_fname = sprintf('%s\\data\\ISIMIP\\Modelregionallocation\\Country_to_REMINDregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\REMIND\\Solar\\Solar_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellCSP{i},CSP_COE_PL{i},GCMID{i},RCPID{i},TIMEID{i},'CSP')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPV{i},PV_COE_PLh{i},GCMID{i},RCPID{i},TIMEID{i},'PV')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPVres{i},PVres_COE_LF{i},GCMID{i},RCPID{i},TIMEID{i},'PVres')
    
    % TIAM
    C2R_fname = sprintf('%s\\data\\ISIMIP\\Modelregionallocation\\Country_to_TIAMregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\TIAM\\Solar\\Solar_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellCSP{i},CSP_COE_PL{i},GCMID{i},RCPID{i},TIMEID{i},'CSP')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPV{i},PV_COE_PLh{i},GCMID{i},RCPID{i},TIMEID{i},'PV')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPVres{i},PVres_COE_LF{i},GCMID{i},RCPID{i},TIMEID{i},'PVres')
    
    % MESSAGE
    C2R_fname = sprintf('%s\\data\\ISIMIP\\Modelregionallocation\\Country_to_MESSAGEregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\MESSAGE\\Solar\\Solar_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellCSP{i},CSP_COE_PL{i},GCMID{i},RCPID{i},TIMEID{i},'CSP')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPV{i},PV_COE_PLh{i},GCMID{i},RCPID{i},TIMEID{i},'PV')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPVres{i},PVres_COE_LF{i},GCMID{i},RCPID{i},TIMEID{i},'PVres')
    
    % COFFEE
    C2R_fname = sprintf('%s\\data\\ISIMIP\\Modelregionallocation\\Country_to_COFFEEregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\COFFEE\\Solar\\Solar_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellCSP{i},CSP_COE_PL{i},GCMID{i},RCPID{i},TIMEID{i},'CSP')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPV{i},PV_COE_PLh{i},GCMID{i},RCPID{i},TIMEID{i},'PV')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPVres{i},PVres_COE_LF{i},GCMID{i},RCPID{i},TIMEID{i},'PVres')
    
    % TIAM31
    C2R_fname = sprintf('%s\\data\\ISIMIP\\Modelregionallocation\\Country_to_TIAM31Rregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\TIAM31R\\Solar\\Solar_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellCSP{i},CSP_COE_PL{i},GCMID{i},RCPID{i},TIMEID{i},'CSP')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPV{i},PV_COE_PLh{i},GCMID{i},RCPID{i},TIMEID{i},'PV')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPVres{i},PVres_COE_LF{i},GCMID{i},RCPID{i},TIMEID{i},'PVres')
    
%     % WITCH
%     C2R_fname = sprintf('%s\\data\\ISIMIP\\Modelregionallocation\\Country_to_WITCHRregion_IsoCode.csv',root);
%     scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\WITCH\\Solar\\Solar_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
%     
%     modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellCSP{i},CSP_COE_PL{i},GCMID{i},RCPID{i},TIMEID{i},'CSP')
%     modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPV{i},PV_COE_PLh{i},GCMID{i},RCPID{i},TIMEID{i},'PV')
%     modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPVres{i},PVres_COE_LF{i},GCMID{i},RCPID{i},TIMEID{i},'PVres')
    
    % TIAM country division
    C2R_fname = sprintf('%s\\data\\ISIMIP\\Modelregionallocation\\GISO_country_division.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\TIAM_country\\Solar\\Solar_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellCSP{i},CSP_COE_PL{i},GCMID{i},RCPID{i},TIMEID{i},'CSP')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPV{i},PV_COE_PLh{i},GCMID{i},RCPID{i},TIMEID{i},'PV')
    modelregoutput(root,C2R_fname,scenlib_IAM,AnnTechPotCellPVres{i},PVres_COE_LF{i},GCMID{i},RCPID{i},TIMEID{i},'PVres')
    
    %% step 6: Deal with IMAGE output
    disp('Output')
    
    scenlib = sprintf('Y:\\Kennisbasis\\IMAGE\\model\\Users\\David\\Pojects\\ISIMIP2E\\ISIMIP2E\\2_TIMER\\scenlib\\TIMER_2015\\ISIMIP2E\\maps\\Solar\\Solar_%s_%s_%s',GCMID{i},RCPID{i},TIMEID{i});
    
    matpath = fullfile(scenlib, sprintf(''));
    if ~isdir(matpath)
        mkdir(matpath);
    end
%     matpath = fullfile(scenlib, sprintf('\\maps'));
%     if ~isdir(matpath)
%         mkdir(matpath);
%     end
    
    pathname = fileparts(scenlib);
    
    % CostCurves
    
    % PV
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\curves\\CostCurveSmthPV.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.PVCostCurveSmth[27](oo) = [');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,CostCurveSmthPV{i},'-append','delimiter','\t');
%     dlmwrite(file,']','-append');
%     dlmwrite(file,';','-append');
%     
%     % CSP
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\curves\\CostCurveSmthCSP.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.CSPCostCurveSmth[27](oo) = [');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,CostCurveSmthCSP{i},'-append','delimiter','\t');
%     dlmwrite(file,']','-append');
%     dlmwrite(file,';','-append');
%     
%     % PVres
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\curves\\CostCurveSmthPVres.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.PVresCostCurveSmth[27](oo) = [');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,CostCurveSmthPVres{i},'-append','delimiter','\t');
%     dlmwrite(file,']','-append');
%     dlmwrite(file,';','-append');
%     
%     % MaxProd
%     
%     % PV
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\curves\\MaxProdPV.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.MaxProdPV[27] = ');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,MaxProdPV{i},'-append','delimiter','\t');
%     dlmwrite(file,';','-append');
%     
%     % CSP
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\curves\\MaxProdCSP.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.MaxProdCSP[27] = ');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,MaxProdCSP{i},'-append','delimiter','\t');
%     dlmwrite(file,';','-append');
%     
%     % PVres
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\curves\\MaxProdPVres.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.MaxProdPVres[27] = ');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,MaxProdPVres{i},'-append','delimiter','\t');
%     dlmwrite(file,';','-append');
%     
%     % Load factor decline curves
%     
%     % PV
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\curves\\LoadCurveSmthPV.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.LFCurveSmthPV[27](oo) = [');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,LFCurveSmthPV{i},'-append','delimiter','\t');
%     dlmwrite(file,']','-append');
%     dlmwrite(file,';','-append');
%     
%     % CSP
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\curves\\LoadCurveSmthCSP.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.LFCurveSmthCSP[27](oo) = [');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,LFCurveSmthCSP{i},'-append','delimiter','\t');
%     dlmwrite(file,']','-append');
%     dlmwrite(file,';','-append');
%     
%     % PVres
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\curves\\LoadCurveSmthPVres.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.LFCurveSmthPVres[27](oo) = [');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,LFCurveSmthPVres{i},'-append','delimiter','\t');
%     dlmwrite(file,']','-append');
%     dlmwrite(file,';','-append');
    
    % Maps
    txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.2f\nNODATA_value\t%d',720,360,-180,-90,0.5,-99);
    
    % PV theoretical
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\TheoPotPV.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    AnnPotCellPV{i}(~imagemask)=-99;
    dlmwrite(file,AnnPotCellPV{i},'-append','delimiter',' '); %kWh / cell / y
    %figure(1);clf;imagesc(AnnPotCellPV{i});axis image
    
    % PV geographic
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\GeoPotPV.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    PV_GeoPot_Cell{i}{13}(~imagemask)=-99;
    dlmwrite(file,PV_GeoPot_Cell{i}{13},'-append','delimiter',' '); %kWh / cell / y
    %figure(2);clf;imagesc(PV_GeoPot_Cell{i}{13});axis image
    
    % PV technical
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\TechPotPV.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    AnnTechPotCellPV{i}(~imagemask)=-99;
    dlmwrite(file,AnnTechPotCellPV{i},'-append','delimiter',' '); % kWh / cell / y
    %figure(3);clf;imagesc(AnnTechPotCellPV{i});axis image

%     if i==6 %for show PhD mid chapter hadgem techpot
%         % PV technical
%         file = fullfile(pathname, sprintf('\\xPot_Mean_change_maps\\Techpot_PhD_midchapter\\data\\Solar_%s_%s_%s_TechPotPV.asc',GCMID{i},RCPID{i},TIMEID{i}));
%         dt=AnnTechPotCellPV{i}*3.6e-9;
%         fprintf('PV techpot %0.2f EJ\n',sum(dt(dt~=-99))*1e-3)
%         dt(dt>90)=90;
%         dt(~imagemask)=-99;
%         dlmwrite(file,txt,'');
%         dlmwrite(file,dt,'-append','delimiter',' '); % PJ / cell / y
%         %figure(2);clf;imagesc(dt);axis image;colorbar
%         
%         % PV economic technical
%         file = fullfile(pathname, sprintf('\\xPot_Mean_change_maps\\Techpot_PhD_midchapter\\data\\Solar_%s_%s_%s_EconPotPV.asc',GCMID{i},RCPID{i},TIMEID{i}));
%         dt=AnnTechPotCellPV{i}*3.6e-9;
%         dt(PV_COE_LF{i}>0.1)=0;
%         fprintf('PV econpot %0.2f EJ\n',sum(dt(dt~=-99))*1e-3)
%         dt(dt>90)=90;
%         dt(~imagemask)=-99;
%         dlmwrite(file,txt,'');
%         dlmwrite(file,dt,'-append','delimiter',' '); % PJ / cell / y
%         %figure(2);clf;imagesc(dt);axis image;colorbar
%     end

    % PV COE
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\COEPV.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    PV_COE_PLh{i}(~imagemask)=-99;
    PV_COE_PLh{i}(find(isnan(PV_COE_PLh{i}(:))))=-99;
    PV_COE_PLh{i}(find(isinf(PV_COE_PLh{i}(:))))=-99;
    PV_COE_PLh{i}(find(PV_COE_PLh{i}>0.5))=0.5;
    dlmwrite(file,PV_COE_PLh{i},'-append','delimiter',' '); % $2010/kWh
    %PV_COE_PLh{i}(PV_COE_PLh{i}==-99)=NaN;
    %figure(4);clf;imagesc(PV_COE_PLh{i});axis image
    
    % CSP theoretical
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\TheoPotCSP.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    AnnPotCellDNI{i}(~imagemask)=-99;
    dlmwrite(file,AnnPotCellDNI{i},'-append','delimiter',' '); %kWh / cell / y
    %figure(1);clf;imagesc(AnnPotCellDNI{i});axis image
    
    % CSP geographic
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\GeoPotCSP.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    CSP_GeoPot_Cell{i}{13}(~imagemask)=-99;
    dlmwrite(file,CSP_GeoPot_Cell{i}{13},'-append','delimiter',' '); %kWh / cell / y
    %figure(1);clf;imagesc(CSP_GeoPot_Cell{i}{13});axis image
    
    % CSP technical
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\TechPotCSP.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    AnnTechPotCellCSP{i}(~imagemask)=-99;
    dlmwrite(file,AnnTechPotCellCSP{i},'-append','delimiter',' ');% kWh / cell / y
    %figure(1);clf;imagesc(AnnTechPotCellCSP{i});axis image
    
	%%
%     if i==6 %for show PhD mid chapter hadgem techpot
%         % CSP technical
%         file = fullfile(pathname, sprintf('\\xPot_Mean_change_maps\\Techpot_PhD_midchapter\\data\\Solar_%s_%s_%s_TechPotCSP.asc',GCMID{i},RCPID{i},TIMEID{i}));
%         dt=AnnTechPotCellCSP{i}*3.6e-9;
%         fprintf('CSP techpot %0.2f EJ\n',sum(dt(dt~=-99))*1e-3)
%         dt(dt>90)=90;
%         dt(~imagemask)=-99;
%         dlmwrite(file,txt,'');
%         dlmwrite(file,dt,'-append','delimiter',' '); % PJ / cell / y
%         %figure(2);clf;imagesc(dt);axis image;colorbar
%         
%         
%         % CSP economic technical
%         file = fullfile(pathname, sprintf('\\xPot_Mean_change_maps\\Techpot_PhD_midchapter\\data\\Solar_%s_%s_%s_EconPotCSP.asc',GCMID{i},RCPID{i},TIMEID{i}));
%         dt=AnnTechPotCellCSP{i}*3.6e-9;
%         dt(CSP_COE_LF{i}>0.1)=0;
%         fprintf('CSP econpot %0.2f EJ\n',sum(dt(dt~=-99))*1e-3)
%         dt(dt>90)=90;
%         dt(~imagemask)=-99;
%         dlmwrite(file,txt,'');
%         dlmwrite(file,dt,'-append','delimiter',' '); % PJ / cell / y
%         %figure(2);clf;imagesc(dt);axis image;colorbar
%     end
    %%
    % CSP COE
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\COECSP.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    CSP_COE_PL{i}(~imagemask)=-99;
    CSP_COE_PL{i}(find(isnan(CSP_COE_PL{i}(:))))=99;
    CSP_COE_PL{i}(find(isinf(CSP_COE_PL{i}(:))))=99;
    CSP_COE_PL{i}(find(CSP_COE_PL{i}>0.5))=0.5;
    dlmwrite(file,CSP_COE_PL{i},'-append','delimiter',' ');  % $2010/kWh
    %CSP_COE_PL{i}(CSP_COE_PL{i}==-99)=NaN;
    %figure(1);clf;imagesc(CSP_COE_PL{i});axis image
    
    % PVres theoretical
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\TheoPotPVres.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    PVres_TheoPot{i}{13}(~imagemask)=-99;
    dlmwrite(file,PVres_TheoPot{i}{13},'-append','delimiter',' '); % kWh / cell / y
    %figure(1);clf;imagesc(PVres_TheoPot{i}{13});axis image
    
    % PVres geographic
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\GeoPotPVres.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    PVres_GeoPot{i}{13}(~imagemask)=-99;
    dlmwrite(file,PVres_GeoPot{i}{13},'-append','delimiter',' '); % kWh / cell / y
    %figure(1);clf;imagesc(PVres_GeoPot{i}{13});axis image
    
    % PVres technical
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\TechPotPVres.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    AnnTechPotCellPVres{i}(~imagemask)=-99;
    dlmwrite(file,AnnTechPotCellPVres{i},'-append','delimiter',' ');% kWh / cell / y
    %figure(1);clf;imagesc(AnnTechPotCellPVres{i});axis image
    
    % PVres COE
    file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\COEPVres.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    PVres_COE_LF{i}(~imagemask)=-99;
    PVres_COE_LF{i}(find(isnan(PVres_COE_LF{i}(:))))=99;
    PVres_COE_LF{i}(find(isinf(PVres_COE_LF{i}(:))))=99;
    PVres_COE_LF{i}(find(PVres_COE_LF{i}>0.5))=0.5;
    dlmwrite(file,PVres_COE_LF{i},'-append','delimiter',' '); % $2010/kWh
    %figure(1);clf;imagesc(PVres_COE_LF{i});axis image
    
%     %% input sce-file
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\scenario_input.sce',GCMID{i},RCPID{i},TIMEID{i}));
%     txt1 = sprintf('DIRECTORY("../scenlib/$1/curves");');
%     txt2 = sprintf('FILE("CostCurveSmthPV.dat","r")     =  main.em.ep.isod.CostCurveSmthPV;');
%     txt3 = sprintf('FILE("CostCurveSmthCSP.dat","r")    =  main.em.ep.isod.CostCurveSmthCSP;');
%     txt4 = sprintf('FILE("CostCurveSmthPVres.dat","r")  =  main.em.ep.isod.CostCurveSmthPVres;');
%     txt5 = sprintf('FILE("MaxProdPV.dat","r")           =  main.em.ep.isod.MaxProdPV;');
%     txt6 = sprintf('FILE("MaxProdCSP.dat","r")          =  main.em.ep.isod.MaxProdCSP;');
%     txt7 = sprintf('FILE("MaxProdPVres.dat","r")        =  main.em.ep.isod.MaxProdPVres;');
%     txt8 = sprintf('FILE("LoadCurveSmthPV.dat","r")     =  main.em.ep.isod.LoadCurveSmthPV;');
%     txt9 = sprintf('FILE("LoadCurveSmthCSP.dat","r")    =  main.em.ep.isod.LoadCurveSmthCSP;');
%     txt10 = sprintf('FILE("LoadCurveSmthPVres.dat","r")  =  main.em.ep.isod.LoadCurveSmthPVres;');
%     
%     dlmwrite(file,txt1,'')
%     dlmwrite(file,txt2,'-append','delimiter','')
%     dlmwrite(file,txt3,'-append','delimiter','')
%     dlmwrite(file,txt4,'-append','delimiter','')
%     dlmwrite(file,txt5,'-append','delimiter','')
%     dlmwrite(file,txt6,'-append','delimiter','')
%     dlmwrite(file,txt7,'-append','delimiter','')
%     dlmwrite(file,txt8,'-append','delimiter','')
%     dlmwrite(file,txt9,'-append','delimiter','')
%     dlmwrite(file,txt10,'-append','delimiter','')
%     
%     %% output sce-file
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\scenario_output.sce',GCMID{i},RCPID{i},TIMEID{i}));
%     txt1 = sprintf('DIRECTORY("../outputlib/$1");');
%     txt2 = sprintf('FILE("ElecProdSpec.out", "w")			= main.em.ep.opr.ElecProdSpec2;');
%     txt3 = sprintf('FILE("GCap.out","w")                            = main.em.ep.gcp.GCap;');
%     txt4 = sprintf('FILE("CO2Spec.out","w") 			= main.em.CO2Spec;');
%     txt5 = sprintf('FILE("tpesEXT.out","w")                         = main.em.TPEStotalEXT;');
%     txt6 = sprintf('FILE("TotalCostPerkWhNew.out", "w")		= main.em.ep.epc.TotalCostPerkWhNew;');
%     txt7 = sprintf('FILE("TotalCostPerkWhAvg.out", "w")		= main.em.ep.epc.TotalCostPerkWhAvg;');
%     txt8 = sprintf('FILE("StorLossISOTot.out", "w")			= main.em.ep.isp.StorLossISOTot;');
%     txt9 = sprintf('FILE("StorLossFracISO.out", "w")		= main.em.ep.isp.StorLossFracISO;');
%     txt10 = sprintf('FILE("ISODeple.out", "w")                       = main.em.ep.isod.ISODeple;');
%     txt11 = sprintf('FILE("ISOAddOnNew.out", "w")			= main.em.ep.epc.ISOAddOnNew;');
%     txt12 = sprintf('FILE("ISOAddOnAvg.out", "w")			= main.em.ep.epc.ISOAddOnAvg;');
%     txt13 = sprintf('FILE("TotalCostperkWhAvgPVres.out","w")         = main.em.ep.epc.TotalCostperkWhAvgPVres;');
%     txt14 = sprintf('FILE("BLF_Prod.out","w") 			= main.em.bio.BCap.BLFProdTot;');
%     txt15 = sprintf('FILE("BSF_Prod.out","w") 			= main.em.bio.BCap.BSFProdTot;');
%     txt16 = sprintf('FILE("tpes.out","w")                         	= main.em.TPEStotal;');
%     
%     dlmwrite(file,txt1,'')
%     dlmwrite(file,txt2,'-append','delimiter','')
%     dlmwrite(file,txt3,'-append','delimiter','')
%     dlmwrite(file,txt4,'-append','delimiter','')
%     dlmwrite(file,txt5,'-append','delimiter','')
%     dlmwrite(file,txt6,'-append','delimiter','')
%     dlmwrite(file,txt7,'-append','delimiter','')
%     dlmwrite(file,txt8,'-append','delimiter','')
%     dlmwrite(file,txt9,'-append','delimiter','')
%     dlmwrite(file,txt10,'-append','delimiter','')
%     dlmwrite(file,txt11,'-append','delimiter','')
%     dlmwrite(file,txt12,'-append','delimiter','')
%     dlmwrite(file,txt13,'-append','delimiter','')
%     dlmwrite(file,txt14,'-append','delimiter','')
%     dlmwrite(file,txt15,'-append','delimiter','')
%     dlmwrite(file,txt16,'-append','delimiter','')
%     
%     %% settings bat-file
%     file = fullfile(pathname, sprintf('\\Solar_%s_%s_%s\\scenario_settings.bat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt1 = sprintf('REM Scenarios settings');
%     dlmwrite(file,txt1,'')
    
end

%% Time dependent output for TIMER

tv=[1971 2000 2050 2085 2100];
% for time depended curves
jv=[1 1 2 3 3 1 1 4 5 5 ... % GFDL-ESM2M RCP26 hist 2050 2085 2100; RCP60 hist 2050 2085 2100
    6 6 7 8 8 6 6 9 10 10 ... % HadGEM-ES RCP26 hist 2050 2085 2100; RCP60 hist 2050 2085 2100
    11 11 12 13 13 11 11 14 15 15 ... %IPSL-CM5A-LR RCP26 hist 2050 2085 2100; RCP60 hist 2050 2085 2100
    16 16 17 18 18 16 16 19 20 20]; % MIROC5 RCP26 hist 2050 2085 2100; RCP60 hist 2050 2085 2100
j=0;
k=0;
l=0;
%Prep order curves
for i=1:numel(jv)
    j=j+1;
    if j>5; j=1; end
    if j==1; k=k+1; l=0; end
    %fprintf('t=%d j=%d jv=%d k=%d l=%d\n',i,j,jv(i),k,l)
    
    for c=1:101 %cost-curve steps
        for R=1:27 %region
            l=l+1;
            if c==1 && R==1
                CostCurveSmthPVt{k}(l)=tv(j);
                CostCurveSmthCSPt{k}(l)=tv(j);
                CostCurveSmthPVrest{k}(l)=tv(j);
                LFCurveSmthPVt{k}(l)=tv(j);
                LFCurveSmthCSPt{k}(l)=tv(j);
                LFCurveSmthPVrest{k}(l)=tv(j);
                l=l+1;
                CostCurveSmthPVt{k}(l)=CostCurveSmthPV{jv(i)}(c,R+1);
                CostCurveSmthCSPt{k}(l)=CostCurveSmthCSP{jv(i)}(c,R+1);
                CostCurveSmthPVrest{k}(l)=CostCurveSmthPVres{jv(i)}(c,R+1);
                LFCurveSmthPVt{k}(l)=LFCurveSmthPV{jv(i)}(c,R+1);
                LFCurveSmthCSPt{k}(l)=LFCurveSmthCSP{jv(i)}(c,R+1);
                LFCurveSmthPVrest{k}(l)=LFCurveSmthPVres{jv(i)}(c,R+1);
            else
                CostCurveSmthPVt{k}(l)=CostCurveSmthPV{jv(i)}(c,R+1);
                CostCurveSmthCSPt{k}(l)=CostCurveSmthCSP{jv(i)}(c,R+1);
                CostCurveSmthPVrest{k}(l)=CostCurveSmthPVres{jv(i)}(c,R+1);
                LFCurveSmthPVt{k}(l)=LFCurveSmthPV{jv(i)}(c,R+1);
                LFCurveSmthCSPt{k}(l)=LFCurveSmthCSP{jv(i)}(c,R+1);
                LFCurveSmthPVrest{k}(l)=LFCurveSmthPVres{jv(i)}(c,R+1);
            end
        end
    end
end

% prep order maxprod
j=0;
k=0;
l=0;
for i=1:numel(jv)
    j=j+1;
    if j>5; j=1; end
    if j==1; k=k+1; l=0; end
    %fprintf('t=%d j=%d jv=%d k=%d l=%d\n',i,j,jv(i),k,l)
    
    for R=1:27 %region
        l=l+1;
        if R==1
            MaxProdPVt{k}(l)=tv(j);
            MaxProdCSPt{k}(l)=tv(j);
            MaxProdPVrest{k}(l)=tv(j);
            l=l+1;
            MaxProdPVt{k}(l)=MaxProdPV{jv(i)}(R);
            MaxProdCSPt{k}(l)=MaxProdCSP{jv(i)}(R);
            MaxProdPVrest{k}(l)=MaxProdPVres{jv(i)}(R);
        else
            MaxProdPVt{k}(l)=MaxProdPV{jv(i)}(R);
            MaxProdCSPt{k}(l)=MaxProdCSP{jv(i)}(R);
            MaxProdPVrest{k}(l)=MaxProdPVres{jv(i)}(R);
        end
    end
end


% Writing
GCMIDt = {'GFDL-ESM2M' 'HADGEM2-ES' 'IPSL-CM5A-LR' 'MIROC5'};
RCPIDt = {'RCP26' 'RCP60'};
c=0;
for i=1:4
    for j=1:2
        c=c+1;
        scenlib = sprintf('Y:\\Kennisbasis\\IMAGE\\model\\users\\david\\Pojects\\ISIMIP2E\\ISIMIP2E\\2_TIMER\\scenlib\\TIMER_2015\\ISIMIP2E\\Solar_%s_%s',GCMIDt{i},RCPIDt{j});
        matpath = fullfile(scenlib, sprintf('\\curves'));
        if ~isdir(matpath); mkdir(matpath); end
        pathname = fileparts(scenlib);
        
        % PV cost
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\CostCurveSmthPV.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.PVCostCurveSmth[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,CostCurveSmthPVt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % CSP cost
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\CostCurveSmthCSP.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.CSPCostCurveSmth[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,CostCurveSmthCSPt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % PVres cost
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\CostCurveSmthPVres.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.PVresCostCurveSmth[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,CostCurveSmthPVrest{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % PV MaxProd
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\MaxProdPV.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.MaxProdPV[27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,MaxProdPVt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % CSP MaxProd
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\MaxProdCSP.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.MaxProdCSP[27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,MaxProdCSPt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % PVres MaxProd
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\MaxProdPVres.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.MaxProdPVres[27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,MaxProdPVrest{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % PV load
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\LoadCurveSmthPV.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.LFCurveSmthPV[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,LFCurveSmthPVt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % CSP load
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\LoadCurveSmthCSP.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.LFCurveSmthCSP[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,LFCurveSmthCSPt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % PVres load
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\LoadCurveSmthPVres.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.LFCurveSmthPVres[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,LFCurveSmthPVrest{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        %% input sce-file
        % Keep format. in the sce-file the format is correct
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\scenario_input.sce',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('DIRECTORY("../scenlib/$1/curves");');
        txt2 = sprintf('FILE("CostCurveSmthPV.dat","r")	=  main.em.ep.isod.CostCurveSmthPV;');
        txt3 = sprintf('FILE("CostCurveSmthCSP.dat","r")	=  main.em.ep.isod.CostCurveSmthCSP;');
        txt3 = sprintf('FILE("CostCurveSmthPVres.dat","r")	=  main.em.ep.isod.CostCurveSmthPVres;');
        txt4 = sprintf('FILE("MaxProdPV.dat","r")              =  main.em.ep.isod.MaxProdPV;');
        txt5 = sprintf('FILE("MaxProdCSP.dat","r")             =  main.em.ep.isod.MaxProdCSP;');
        txt5 = sprintf('FILE("MaxProdPVres.dat","r")             =  main.em.ep.isod.MaxProdPVres;');
        txt6 = sprintf('FILE("LoadCurveSmthPV.dat","r")	=  main.em.ep.isod.LoadCurveSmthPV;');
        txt7 = sprintf('FILE("LoadCurveSmthCSP.dat","r")	=  main.em.ep.isod.LoadCurveSmthCSP;');
        txt7 = sprintf('FILE("LoadCurveSmthPVres.dat","r")	=  main.em.ep.isod.LoadCurveSmthPVres;');
        
        dlmwrite(file,txt1,'')
        dlmwrite(file,txt2,'-append','delimiter','')
        dlmwrite(file,txt3,'-append','delimiter','')
        dlmwrite(file,txt4,'-append','delimiter','')
        dlmwrite(file,txt5,'-append','delimiter','')
        dlmwrite(file,txt6,'-append','delimiter','')
        dlmwrite(file,txt7,'-append','delimiter','')
        
        %% output sce-file
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\scenario_output.sce',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('DIRECTORY("../outputlib/$1");');
        txt2 = sprintf('FILE("ElecProdSpec.out", "w")			= main.em.ep.opr.ElecProdSpec2;');
        txt3 = sprintf('FILE("GCap.out","w")                            = main.em.ep.gcp.GCap;');
        txt4 = sprintf('FILE("CO2Spec.out","w") 			= main.em.CO2Spec;');
        txt5 = sprintf('FILE("tpesEXT.out","w")                         = main.em.TPEStotalEXT;');
        txt6 = sprintf('FILE("TotalCostPerkWhNew.out", "w")		= main.em.ep.epc.TotalCostPerkWhNew;');
        txt7 = sprintf('FILE("TotalCostPerkWhAvg.out", "w")		= main.em.ep.epc.TotalCostPerkWhAvg;');
        txt8 = sprintf('FILE("StorLossISOTot.out", "w")			= main.em.ep.isp.StorLossISOTot;');
        txt9 = sprintf('FILE("StorLossFracISO.out", "w")		= main.em.ep.isp.StorLossFracISO;');
        txt10 = sprintf('FILE("ISODeple.out", "w")                       = main.em.ep.isod.ISODeple;');
        txt11 = sprintf('FILE("ISOAddOnNew.out", "w")			= main.em.ep.epc.ISOAddOnNew;');
        txt12 = sprintf('FILE("ISOAddOnAvg.out", "w")			= main.em.ep.epc.ISOAddOnAvg;');
        txt13 = sprintf('FILE("TotalCostperkWhAvgPVres.out","w")         = main.em.ep.epc.TotalCostperkWhAvgPVres;');
        txt14 = sprintf('FILE("BLF_Prod.out","w") 			= main.em.bio.BCap.BLFProdTot;');
        txt15 = sprintf('FILE("BSF_Prod.out","w") 			= main.em.bio.BCap.BSFProdTot;');
        txt16 = sprintf('FILE("tpes.out","w")                         	= main.em.TPEStotal;');
        
        dlmwrite(file,txt1,'')
        dlmwrite(file,txt2,'-append','delimiter','')
        dlmwrite(file,txt3,'-append','delimiter','')
        dlmwrite(file,txt4,'-append','delimiter','')
        dlmwrite(file,txt5,'-append','delimiter','')
        dlmwrite(file,txt6,'-append','delimiter','')
        dlmwrite(file,txt7,'-append','delimiter','')
        dlmwrite(file,txt8,'-append','delimiter','')
        dlmwrite(file,txt9,'-append','delimiter','')
        dlmwrite(file,txt10,'-append','delimiter','')
        dlmwrite(file,txt11,'-append','delimiter','')
        dlmwrite(file,txt12,'-append','delimiter','')
        dlmwrite(file,txt13,'-append','delimiter','')
        dlmwrite(file,txt14,'-append','delimiter','')
        dlmwrite(file,txt15,'-append','delimiter','')
        dlmwrite(file,txt16,'-append','delimiter','')
        
        %% settings bat-file
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\scenario_settings.bat',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('REM Scenarios settings');
        dlmwrite(file,txt1,'')
        
    end
end

%% Time dependent output for TIMER without climate change impacts, just historic values
clear tv jv CostCurveSmthPVt CostCurveSmthCSPt CostCurveSmthPVrest LFCurveSmthPVt LFCurveSmthCSPt LFCurveSmthPVrest
clear MaxProdPVt MaxProdCSPt MaxProdPVrest
tv=[1971 2000 2050 2100];
% for historic only curves
jv=[1 1 1 1 1 1 1 1 ... % GFDL-ESM2M RCP26 hist hist hist; RCP60 hist hist hist
    6 6 6 6 6 6 6 6 ... % HadGEM-ES RCP26 hist hist hist; RCP60 hist hist hist
    11 11 11 11 11 11 11 11 ... %IPSL-CM5A-LR RCP26 hist hist hist; RCP60 hist hist hist
    16 16 16 16 16 16 16 16]; % MIROC5 RCP26 hist hist hist; RCP60 hist hist hist
j=0;
k=0;
l=0;
%Prep order curves
for i=1:numel(jv)
    j=j+1;
    if j>4; j=1; end
    if j==1; k=k+1; l=0; end
    %fprintf('t=%d j=%d jv=%d k=%d l=%d\n',i,j,jv(i),k,l)
    
    for c=1:101 %cost-curve steps
        for R=1:27 %region
            l=l+1;
            if c==1 && R==1
                CostCurveSmthPVt{k}(l)=tv(j);
                CostCurveSmthCSPt{k}(l)=tv(j);
                CostCurveSmthPVrest{k}(l)=tv(j);
                LFCurveSmthPVt{k}(l)=tv(j);
                LFCurveSmthCSPt{k}(l)=tv(j);
                LFCurveSmthPVrest{k}(l)=tv(j);
                l=l+1;
                CostCurveSmthPVt{k}(l)=CostCurveSmthPV{jv(i)}(c,R+1);
                CostCurveSmthCSPt{k}(l)=CostCurveSmthCSP{jv(i)}(c,R+1);
                CostCurveSmthPVrest{k}(l)=CostCurveSmthPVres{jv(i)}(c,R+1);
                LFCurveSmthPVt{k}(l)=LFCurveSmthPV{jv(i)}(c,R+1);
                LFCurveSmthCSPt{k}(l)=LFCurveSmthCSP{jv(i)}(c,R+1);
                LFCurveSmthPVrest{k}(l)=LFCurveSmthPVres{jv(i)}(c,R+1);
            else
                CostCurveSmthPVt{k}(l)=CostCurveSmthPV{jv(i)}(c,R+1);
                CostCurveSmthCSPt{k}(l)=CostCurveSmthCSP{jv(i)}(c,R+1);
                CostCurveSmthPVrest{k}(l)=CostCurveSmthPVres{jv(i)}(c,R+1);
                LFCurveSmthPVt{k}(l)=LFCurveSmthPV{jv(i)}(c,R+1);
                LFCurveSmthCSPt{k}(l)=LFCurveSmthCSP{jv(i)}(c,R+1);
                LFCurveSmthPVrest{k}(l)=LFCurveSmthPVres{jv(i)}(c,R+1);
            end
            
        end
        
    end
end

% prep order maxprod
j=0;
k=0;
l=0;
for i=1:numel(jv)
    j=j+1;
    if j>4; j=1; end
    if j==1; k=k+1; l=0; end
    %fprintf('t=%d j=%d jv=%d k=%d l=%d\n',i,j,jv(i),k,l)
    
    for R=1:27 %region
        l=l+1;
        if R==1
            MaxProdPVt{k}(l)=tv(j);
            MaxProdCSPt{k}(l)=tv(j);
            MaxProdPVrest{k}(l)=tv(j);
            l=l+1;
            MaxProdPVt{k}(l)=MaxProdPV{jv(i)}(R);
            MaxProdCSPt{k}(l)=MaxProdCSP{jv(i)}(R);
            MaxProdPVrest{k}(l)=MaxProdPVres{jv(i)}(R);
        else
            MaxProdPVt{k}(l)=MaxProdPV{jv(i)}(R);
            MaxProdCSPt{k}(l)=MaxProdCSP{jv(i)}(R);
            MaxProdPVrest{k}(l)=MaxProdPVres{jv(i)}(R);
        end
    end
end


% Writing
GCMIDt = {'GFDL-ESM2M' 'HADGEM2-ES' 'IPSL-CM5A-LR' 'MIROC5'};
RCPIDt = {'hist'};
c=0;
for i=1:4
    for j=1
        c=c+2;
        scenlib = sprintf('Y:\\Kennisbasis\\IMAGE\\model\\users\\david\\Pojects\\ISIMIP2E\\ISIMIP2E\\2_TIMER\\scenlib\\TIMER_2015\\ISIMIP2E\\Solar_%s_%s',GCMIDt{i},RCPIDt{j});
        matpath = fullfile(scenlib, sprintf('\\curves'));
        if ~isdir(matpath); mkdir(matpath); end
        pathname = fileparts(scenlib);
        
        % PV cost
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\CostCurveSmthPV.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.PVCostCurveSmth[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,CostCurveSmthPVt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % CSP cost
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\CostCurveSmthCSP.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.CSPCostCurveSmth[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,CostCurveSmthCSPt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % PVres cost
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\CostCurveSmthPVres.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.PVresCostCurveSmth[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,CostCurveSmthPVrest{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % PV MaxProd
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\MaxProdPV.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.MaxProdPV[27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,MaxProdPVt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % CSP MaxProd
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\MaxProdCSP.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.MaxProdCSP[27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,MaxProdCSPt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % PVres MaxProd
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\MaxProdPVres.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.MaxProdPVres[27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,MaxProdPVrest{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % PV load
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\LoadCurveSmthPV.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.LFCurveSmthPV[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,LFCurveSmthPVt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % CSP load
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\LoadCurveSmthCSP.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.LFCurveSmthCSP[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,LFCurveSmthCSPt{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % PVres load
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\curves\\LoadCurveSmthPVres.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.LFCurveSmthPVres[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,LFCurveSmthPVrest{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        %% input sce-file
        % Keep format. in the sce-file the format is correct
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\scenario_input.sce',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('DIRECTORY("../scenlib/$1/curves");');
        txt2 = sprintf('FILE("CostCurveSmthPV.dat","r")	=  main.em.ep.isod.CostCurveSmthPV;');
        txt3 = sprintf('FILE("CostCurveSmthCSP.dat","r")	=  main.em.ep.isod.CostCurveSmthCSP;');
        txt3 = sprintf('FILE("CostCurveSmthPVres.dat","r")	=  main.em.ep.isod.CostCurveSmthPVres;');
        txt4 = sprintf('FILE("MaxProdPV.dat","r")              =  main.em.ep.isod.MaxProdPV;');
        txt5 = sprintf('FILE("MaxProdCSP.dat","r")             =  main.em.ep.isod.MaxProdCSP;');
        txt5 = sprintf('FILE("MaxProdPVres.dat","r")             =  main.em.ep.isod.MaxProdPVres;');
        txt6 = sprintf('FILE("LoadCurveSmthPV.dat","r")	=  main.em.ep.isod.LoadCurveSmthPV;');
        txt7 = sprintf('FILE("LoadCurveSmthCSP.dat","r")	=  main.em.ep.isod.LoadCurveSmthCSP;');
        txt7 = sprintf('FILE("LoadCurveSmthPVres.dat","r")	=  main.em.ep.isod.LoadCurveSmthPVres;');
        
        dlmwrite(file,txt1,'')
        dlmwrite(file,txt2,'-append','delimiter','')
        dlmwrite(file,txt3,'-append','delimiter','')
        dlmwrite(file,txt4,'-append','delimiter','')
        dlmwrite(file,txt5,'-append','delimiter','')
        dlmwrite(file,txt6,'-append','delimiter','')
        dlmwrite(file,txt7,'-append','delimiter','')
        
        %% output sce-file
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\scenario_output.sce',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('DIRECTORY("../outputlib/$1");');
        txt2 = sprintf('FILE("ElecProdSpec.out", "w")			= main.em.ep.opr.ElecProdSpec2;');
        txt3 = sprintf('FILE("GCap.out","w")                            = main.em.ep.gcp.GCap;');
        txt4 = sprintf('FILE("CO2Spec.out","w") 			= main.em.CO2Spec;');
        txt5 = sprintf('FILE("tpesEXT.out","w")                         = main.em.TPEStotalEXT;');
        txt6 = sprintf('FILE("TotalCostPerkWhNew.out", "w")		= main.em.ep.epc.TotalCostPerkWhNew;');
        txt7 = sprintf('FILE("TotalCostPerkWhAvg.out", "w")		= main.em.ep.epc.TotalCostPerkWhAvg;');
        txt8 = sprintf('FILE("StorLossISOTot.out", "w")			= main.em.ep.isp.StorLossISOTot;');
        txt9 = sprintf('FILE("StorLossFracISO.out", "w")		= main.em.ep.isp.StorLossFracISO;');
        txt10 = sprintf('FILE("ISODeple.out", "w")                       = main.em.ep.isod.ISODeple;');
        txt11 = sprintf('FILE("ISOAddOnNew.out", "w")			= main.em.ep.epc.ISOAddOnNew;');
        txt12 = sprintf('FILE("ISOAddOnAvg.out", "w")			= main.em.ep.epc.ISOAddOnAvg;');
        txt13 = sprintf('FILE("TotalCostperkWhAvgPVres.out","w")         = main.em.ep.epc.TotalCostperkWhAvgPVres;');
        txt14 = sprintf('FILE("BLF_Prod.out","w") 			= main.em.bio.BCap.BLFProdTot;');
        txt15 = sprintf('FILE("BSF_Prod.out","w") 			= main.em.bio.BCap.BSFProdTot;');
        txt16 = sprintf('FILE("tpes.out","w")                         	= main.em.TPEStotal;');
        
        dlmwrite(file,txt1,'')
        dlmwrite(file,txt2,'-append','delimiter','')
        dlmwrite(file,txt3,'-append','delimiter','')
        dlmwrite(file,txt4,'-append','delimiter','')
        dlmwrite(file,txt5,'-append','delimiter','')
        dlmwrite(file,txt6,'-append','delimiter','')
        dlmwrite(file,txt7,'-append','delimiter','')
        dlmwrite(file,txt8,'-append','delimiter','')
        dlmwrite(file,txt9,'-append','delimiter','')
        dlmwrite(file,txt10,'-append','delimiter','')
        dlmwrite(file,txt11,'-append','delimiter','')
        dlmwrite(file,txt12,'-append','delimiter','')
        dlmwrite(file,txt13,'-append','delimiter','')
        dlmwrite(file,txt14,'-append','delimiter','')
        dlmwrite(file,txt15,'-append','delimiter','')
        dlmwrite(file,txt16,'-append','delimiter','')
        
        %% settings bat-file
        file = fullfile(pathname, sprintf('\\Solar_%s_%s\\scenario_settings.bat',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('REM Scenarios settings');
        dlmwrite(file,txt1,'')
        
    end
end

%% Visualization

% figure(1);clf;imagesc(CM1_data(:,:,1)'); axis image; colormap(jet)%colormap(parula)
%
% figure(2);clf;imagesc(CM1{13}(:,:)); axis image; colormap(jet)%colormap(parula)
%
% figure(3);clf;
% for m=1:12
%     subplot(3,4,m)
%     imagesc(CM1{m}); axis image; colormap(jet)
% end

