%Data prep
clear all
root = 'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER\CSP_PV_PVres';

%% Image land mask
disp('IMAGE land mask')
fname = sprintf('%s\\data\\imagemask_land.mat', root);
load(fname)

%% NASA DNI
disp('NASA DNI')
fname = sprintf('%s\\data\\CSPnasa_dni_data.dat', root); % kWh / m2 / day
[Avg_DNI_NASA_data] = read_mym(fname);

%% Convert NASA DNI to map
for m=1:13
    Avg_DNI_NASA{m} = zeros(360,720);
end
[nr nc] = size(Avg_DNI_NASA{1});


for m=1:13
    i=0;
    for r=1:nr
        for c=1:nc
            if imagemask(r,c)==0; continue; end;
            i=i+1;
            Avg_DNI_NASA{m}(r,c)=Avg_DNI_NASA_data(i,m); %kwh/m2/day
        end
    end
end

%% Check
% figure(1);clf;
% for m=1:12
%     subplot(3,4,m)
%     imagesc(Avg_DNI_NASA{m});colormap(jet); axis image
% end
%
% for m=1:13; mcell(m) = Avg_DNI_NASA{m}(76,372); end
%
% figure(2);clf; plot(1:13,mcell);


%% NASA horizontal annual
disp('NASA Horizontal')
fname = sprintf('%s\\data\\global_horizontal_radiation_forPV.dat', root); %kWh / m2 / day
[Gl_Horiz_Ann_NASA_data] = read_mym(fname);

%% NASA horizontal data monthly

mname{1} = 'Jan';
mname{2} = 'Feb';
mname{3} = 'Mar';
mname{4} = 'Apr';
mname{5} = 'May';
mname{6} = 'Jun';
mname{7} = 'Jul';
mname{8} = 'Aug';
mname{9} = 'Sep';
mname{10} = 'Oct';
mname{11} = 'Nov';
mname{12} = 'Dec';

for m=1:12
    if m<10
        fname = sprintf('%s\\data\\Gl_Horiz_0%d_%s.asc', root,m, mname{m});
    else
        fname = sprintf('%s\\data\\Gl_Horiz_%d_%s.asc', root,m, mname{m});
    end
    [Gl_Horiz_data{m}] = read_mym(fname);
end

%% Put in map
for m=1:13
    Gl_Horiz_NASA{m} = zeros(360,720); %kwh/m2/day
end

for m=1:13
    i=0;
    for r=1:nr
        for c=1:nc
            if imagemask(r,c)==0; continue; end;
            i=i+1;
            if m<13
                Gl_Horiz_NASA{m}(r,c)=Gl_Horiz_data{m}(i);
            else
                Gl_Horiz_NASA{m}(r,c)=Gl_Horiz_Ann_NASA_data(i);
            end
        end
    end
end

%%
% figure(1);clf;
% for m=1:12
%     subplot(3,4,m)
%     imagesc(Gl_Horiz_NASA{m});colormap(jet); axis image
% end
%
% for m=1:13; mcell(m) = Gl_Horiz_NASA{m}(76,372); end
%
% figure(2);clf; plot(1:13,mcell);

%% IMAGE regions
disp('IMAGE Regions')
fname = sprintf('%s\\data\\region27.dat', root);
[Region_data] = read_mym(fname);

IRegion = zeros(360,720);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        IRegion(r,c)=Region_data(i);
    end
end

% figure(1);clf;imagesc(IRegion);axis image; colormap(jet)

%% Countries
disp('Countries')
fname = sprintf('%s\\data\\gISO2.asc', root);
GISO = zeros(360,720);
GISO(2:end,2:end) = dlmread(fname,',',7,1);
GISO(GISO==-9999)=NaN;
% imagesc(GISO); colorbar

%% Countries
% disp('Countries')
% fname = sprintf('%s\\data\\Countries.mat', root);
% load(fname);
%
% Countries = zeros(360,720);
%
% scalingf = 6; %scaling factor
%
% [nr nc] = size(Countries2);
%
% rb = linspace(scalingf,nr,nr/scalingf);
% cb = linspace(scalingf,nc,nc/scalingf);
%
% Countries = zeros((size(Countries2)/scalingf),'single'); %Weighted Distance map
% [nrW, ncW] = size(Countries);
%
% Countries2=single(Countries2);
%
% r1=1;
% c1=1;
% i=0;
% for r=rb
%     i=i+1;
%     if i==1; r1=1; else r1=rb(i-1)+1; end
%     r2=r;
%     j=0;
%
%     for c=cb
%         j=j+1;
%         if j==1; c1=1; else c1=cb(j-1)+1; end
%         c2=c;
%         Block = Countries2(r1:r2,c1:c2);
%
%         Countries(i,j) = Block(12);
%
%     end
%
% end

% figure(1);clf;imagesc(Countries); axis image;

%% Area regions
disp('Area per cell')
fname = sprintf('%s\\data\\area.dat', root);
[Area_data] = read_mym(fname);

Area = zeros(360,720); %km2
[nr nc] = size(Area);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        Area(r,c)=Area_data(i);
    end
end

% figure(1);clf;imagesc(Area);axis image; colormap(jet); colorbar

%% GLCT new SSPs
disp('GLCT')
fname = sprintf('%s\\data\\GLCT_SSP2.asc', root);
[GLCT_data,t] = read_mym(fname);

% Put in map
GLCT = zeros(360,720); %Landcover 2010
[nr,nc]=size(GLCT);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        GLCT(r,c)=GLCT_data(9,i);
    end
end

% figure(1);clf;imagesc(GLCT);axis image; colormap(prism); colorbar

%% Bioreservce
disp('Bio reserves')
fname = sprintf('%s\\data\\BIORESERVE.dat', root);
[BIOres_data t] = read_mym(fname);

% Put in map
BIOres = zeros(360,720); %Bioreserves [type 0, 1 or 2], value of 0 means available

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        BIOres(r,c)=BIOres_data(i);
    end
end

% figure(1);clf;imagesc(BIOres);axis image; colormap(jet); colorbar

%% Buildup
disp('Build Up')
fname = sprintf('%s\\data\\BUILDUP.dat', root);
[Buildup_data t] = read_mym(fname);

% Put in map
Buildup = zeros(360,720); %Urban areas (Fraction per cell)

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        Buildup(r,c)=Buildup_data(i);
    end
end

% figure(1);clf;imagesc(Buildup);axis image; colormap(jet); colorbar

%% Pop
disp('Population')
fname = sprintf('%s\\data\\totpop.dat', root);
[totpop_data t] = read_mym(fname);
fname = sprintf('%s\\data\\urbpop.dat', root);
[urbpop_data t] = read_mym(fname);
fname = sprintf('%s\\data\\rurpop.dat', root);
[rurpop_data t] = read_mym(fname);

% Put in map
totpop = zeros(360,720);
urbpop = zeros(360,720);
rurpop = zeros(360,720);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        totpop(r,c)=totpop_data(i);
        urbpop(r,c)=urbpop_data(i);
        rurpop(r,c)=rurpop_data(i);
    end
end

% figure(1);clf;imagesc(log(totpop));axis image; colormap(jet); colorbar
% figure(2);clf;imagesc(log(urbpop));axis image; colormap(jet); colorbar
% figure(3);clf;imagesc(log(rurpop));axis image; colormap(jet); colorbar

%% Slope
disp('Slope')
fname = sprintf('%s\\data\\Slope_lessthan_3pcnt.dat', root);
[Slope_data t] = read_mym(fname);

% Put in map
Slope = zeros(360,720); %!Fraction of cell area with slope <3%

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        Slope(r,c)=Slope_data(i);
    end
end

% figure(1);clf;imagesc(Slope);axis image; colormap(jet); colorbar

%% Floorspace and household data
disp('Floorspace and Households')
fname = sprintf('%s\\data\\res_Floorspace.out', root);
[resFloorspace_TURQ t_rFS] = read_mym(fname);

fname = sprintf('%s\\data\\res_Households.out', root);
[NumHHs_TURQ t] = read_mym(fname);

fname = sprintf('%s\\data\\POP_q.out', root);
[Pop_TURQ t] = read_mym(fname);

%% DistanceToLoad
disp('Distance to load')

fname = sprintf('%s\\data\\transmission2010.dat', root);
[DistanceToLoad_data ~] = read_mym(fname);

% Put in map
DistanceToLoad = zeros(360,720);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        if DistanceToLoad_data(1,i)==0; DistanceToLoad(r,c)=1; continue; end;
        DistanceToLoad(r,c)=DistanceToLoad_data(1,i);
    end
end

% DistanceToLoad(DistanceToLoad(:)>1000)=1000;

% figure(1);clf;imagesc(log(totpop));axis image; colormap(jet); colorbar

%% IRRADIANCE from ISIMIP
disp('Solar data prep')
fname = sprintf('%s\\data\\rsds_day_HadGEM2-ES_histor_r1i1p1_EWEMBI_19710101-20001231_mergetime_ymonmean.nc4',root);

% finfo = ncinfo(fname);
% finfo.Variables.Name

CM_data = ncread(fname,'rsds');

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
            CMc{m}(r,c) = CM{m}(r,c) / 41.667; %kwh/m2/day
        end
    end
end

% figure(1);clf;imagesc(CMc{13});colorbar

%% TEMPERATURE FROM ISIMIP
disp('Temperature data prep')
fname = sprintf('%s\\data\\tas_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_19710101-20001231_mergetime_ymonmean.nc4',root);

%finfo = ncinfo(fname);
%finfo.Variables.Name

tas_data = ncread(fname,'tas');

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

%% Wind Speeds FROM ISIMIP
disp('Wind data prep')
fname = sprintf('%s\\data\\sfcWind_day_HadGEM2-ES_histor_r1i1p1_EWEMBI_19710101-20001231_mergetime_ymonmean.nc4',root);
%finfo = ncinfo(fname);
%finfo.Variables.Name

ws_data = ncread(fname,'sfcWind');

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

%% IMAGE regions names
fprintf('IMAGE regions names\n')

IMAGER{1}='Canada';
IMAGER{2}='USA';
IMAGER{3}='Mexico';
IMAGER{4}='Rest Central America';
IMAGER{5}='Brazil';
IMAGER{6}='Rest South America';
IMAGER{7}='Northern Africa';
IMAGER{8}='Western Africa';
IMAGER{9}='Eastern Africa';
IMAGER{10}='Southern Africa';
IMAGER{11}='Western Europe';
IMAGER{12}='Central Europe';
IMAGER{13}='Turkey';
IMAGER{14}='Ukraine +';
IMAGER{15}='Asia-Stan';
IMAGER{16}='Russia +';
IMAGER{17}='Middle East';
IMAGER{18}='India +';
IMAGER{19}='Korea';
IMAGER{20}='China +';
IMAGER{21}='Southeastern Asia';
IMAGER{22}='Indonesia +';
IMAGER{23}='Japan';
IMAGER{24}='Oceania';
IMAGER{25}='Rest S.Asia';
IMAGER{26}='Rest S.Africa';
IMAGER{27}='World';

%% Reading GDPpc ppp
fprintf('GDP pc PPP\n')

fname = sprintf('%s\\data\\GDPpc2010IsoCode.csv', root);
fileID = fopen(fname);
C = textscan(fileID,'%s %s %s %s %s %s %s','Delimiter',';','HeaderLines',1);
fclose(fileID);
for i=1:numel(C{3})
    ISOGDP(i,1) = str2num(C{3}{i})'; %ISO number
    ISOGDP(i,2) = str2num(C{4}{i})'; %GDP pc
    ISOGDP(i,3) = str2num(C{5}{i})'; %GDP MER
    
    CountryNames{i} = C{1}{i};       %Country names
end

%% latmap
fprintf('Latmap\n')

fname = sprintf('%s\\data\\latmap.mat', root);
load(fname)

%% Calculate Betas that convert floor space in roof area
fprintf('Betas floor2roof\n')

fname = sprintf('%s\\data\\floor_space_data.csv', root);
floorspace_data=csvread(fname,1,3); % m2/household for the same image region and the same year as the statistical data
fileID = fopen(fname);
C = textscan(fileID,'%s %s %s %d %d %d','Delimiter',',','HeaderLines',1);
fclose(fileID);
C{1}(:)'; %available countries

fname = sprintf('%s\\data\\floor_data.csv', root); % percentage of households with x amount of floor (1 to 10 floors)
floor_data=csvread(fname,1,1);

fname = sprintf('%s\\data\\nhouseholds_data.csv', root); % percentage of households with x amount of floor (1 to 10 floors)
nhousehold_data=csvread(fname,1,1);

% Calculate total roof area based on statistical data
[nr,nc]=size(floor_data);

for r=1:nr %number of regions
    for f=1:10 %number of floors
        %fprintf('%d\n',f)
        roofarea(r,f) = (floorspace_data(r,3).*floor_data(r,f))/f;
    end
    roofarea(r,11) = sum(roofarea(r,:)); %the average roof area for the region
end

% Calculate which Beta is required to get from floor space to roof area
for r=1:nr %number of regions
    Beta(r) = roofarea(r,11)/floorspace_data(r,3);
end
Beta(nr+1)= mean(Beta);

% Calculate weighted averages for the IMAGE regions with more than one country with statistical data

%Middle East (1=Qatar,7=Jordan,12=Bahrain)
nh_me = nhousehold_data(1,3)+nhousehold_data(7,3)+nhousehold_data(12,3);
Beta_wa(1) = (Beta(1)*nhousehold_data(1,3))/nh_me + (Beta(7)*nhousehold_data(7,3))/nh_me + (Beta(12)*nhousehold_data(12,3))/nh_me;

% Western Europe(8=Finland,9=Greece,10=Switserland,17=Germany,21=Spain)
nh_we = nhousehold_data(8,3)+nhousehold_data(9,3)+nhousehold_data(10,3)+nhousehold_data(17,3)+nhousehold_data(21,3);
Beta_wa(2) = (Beta(8)*nhousehold_data(8,3))/nh_we + (Beta(9)*nhousehold_data(9,3))/nh_we + (Beta(10)*nhousehold_data(10,3))/nh_we ....
    + (Beta(17)*nhousehold_data(17,3))/nh_we + (Beta(21)*nhousehold_data(21,3))/nh_we;

% Eastern Europe(3=Hungary,5=Albania,11=Czech Republic)
nh_ee = nhousehold_data(3,3)+nhousehold_data(5,3)+nhousehold_data(11,3);
Beta_wa(3) = (Beta(3)*nhousehold_data(3,3))/nh_ee + (Beta(5)*nhousehold_data(5,3))/nh_ee + (Beta(11)*nhousehold_data(11,3))/nh_ee;

% World
nh_w = sum(nhousehold_data(:,3));
for r=1:nr
    dt(r) = (Beta(r)*nhousehold_data(r,3))/nh_w;
end
Beta_wa(4) = sum(dt(:));

% Math Beta with IMAGE regions
Beta_fs_ra(1) = Beta(2); %Canada (2=US)
Beta_fs_ra(2) = Beta(2); %US (2=US)
Beta_fs_ra(3) = Beta_wa(4); %Mexico (Beta_wa(4)=weighted average)
Beta_fs_ra(4) = Beta_wa(4); %Rest of Central America (Beta_wa(4)=weighted average)
Beta_fs_ra(5) = Beta_wa(4); %Brazil (Beta_wa(4)=weighted average)
Beta_fs_ra(6) = Beta_wa(4); %Rest of South America (Beta_wa(4)=weighted average)
Beta_fs_ra(7) = Beta_wa(1); % Northern Africa (1=Qatar,7=Jordan,12=Bahrain) 
Beta_fs_ra(8) = Beta(13); % Western Africa (13=Cape Verde)
Beta_fs_ra(9) = Beta(15); % Eastern Africa (15=Mauritius)
Beta_fs_ra(10) = Beta_wa(4); % Southern Africa (Beta_wa(4)=weighted average)
Beta_fs_ra(11) = Beta_wa(2); % Western Europe (8=Finland,9=Greece,10=Switserland,17=Germany,21=Spain)
Beta_fs_ra(12) = Beta_wa(3); % Central Europe (3=Hungary,5=Albania,11=Czech Republic)
Beta_fs_ra(13) = Beta(14); % Turkey (14=Turkey)
Beta_fs_ra(14) = Beta_wa(3); % Ukraine + (3=Hungary,5=Albania,11=Czech Republic)
Beta_fs_ra(15) = Beta(20); % Asia-Stan (20=Kyrgyz Republic)
Beta_fs_ra(16) = Beta(20); % Russia + (20=Kyrgyz Republic)
Beta_fs_ra(17) = Beta_wa(1); % Middle East (1=Qatar,7=Jordan,12=Bahrain)
Beta_fs_ra(18) = Beta(19); % India + (19=Indonesia)
Beta_fs_ra(19) = Beta(6); % Korea (6=Japan)
Beta_fs_ra(20) = Beta_wa(4); % China + (Beta_wa(4)=weighted average)
Beta_fs_ra(21) = Beta(19); % Southeastern Asia (19=Indonesia)
Beta_fs_ra(22) = Beta(19); % Indonesia + (19=Indonesia)
Beta_fs_ra(23) = Beta(6); % Japan (6=Japan)
Beta_fs_ra(24) = Beta(4); % Oceania (4=Australia)
Beta_fs_ra(25) = Beta(19); % Rest S.Asia' (19=Indonesia)
Beta_fs_ra(26) = Beta_wa(4); % Rest S.Africa (Beta_wa(4)=weighted average)
Beta_fs_ra(27) = Beta_wa(4); % World (Beta_wa(4)=weighted average)

%% save
fprintf('Save\n')
% fname = sprintf('%s\\data\\input_data.mat', root);
% load(fname)

matfile = fullfile(root, sprintf('data\\input_data.mat'));
save(matfile,'Avg_DNI_NASA','Gl_Horiz_NASA','IRegion','Area','GLCT','BIOres','Buildup','imagemask','totpop',...
    'urbpop','rurpop','NumHHs_TURQ','resFloorspace_TURQ','t_rFS','Slope','DistanceToLoad','IMAGER',...
    'Pop_TURQ','GISO','ISOGDP','CountryNames','latmap','tasc','wsc','CMc','Beta_fs_ra');

matfile = fullfile(root, sprintf('data\\input_data_ISIMIP.mat'));
save(matfile,'IRegion','Area','GLCT','BIOres','Buildup','imagemask','totpop',...
    'urbpop','rurpop','NumHHs_TURQ','resFloorspace_TURQ','t_rFS','Slope','DistanceToLoad','IMAGER',...
    'Pop_TURQ','GISO','ISOGDP','CountryNames','latmap','Beta_fs_ra');
