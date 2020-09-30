 %% Adjusted DNI
DNI_DailyMin = 3;

for m=1:12
    for r=1:nr
        for c=1:nc
            
            if Avg_DNI_NASA{m}(r,c) < DNI_DailyMin
                Avg_DNI_NASA_adj{m}(r,c)=0;
            else
                Avg_DNI_NASA_adj{m}(r,c)=Avg_DNI_NASA{m}(r,c); %kWh / m2 / day
            end
            
        end
    end
end

%% New annual average
for r=1:nr
    for c=1:nc
        clear Sum_temp
        for m=1:12
            Sum_temp(m) = Avg_DNI_NASA_adj{m}(r,c);
        end
        Avg_DNI_NASA_adj{13}(r,c) = mean(Sum_temp); %kWh / m2 / day
    end
end

%% Insolution kwh/m2/day --> W/m2 conversion

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
    for r=1:nr
        for c=1:nc
            Avg_DNI_NASA_adj_Insol{m}(r,c) = Avg_DNI_NASA_adj{m}(r,c) * 41.667;% 1142; % W / m2 insolation
        end
    end
end

%% DNI Regional and Global Potentials
m2tokm2conv = 1000000;		%1 km2 = 1.000.000 m2
for r=1:nr
    for c=1:nc
        AnnPotCellDNI(r,c) = Avg_DNI_NASA_adj{13}(r,c) * Area(r,c) * m2tokm2conv * 365; %kWh / cell / y
    end
end

% Image reg
for i=1:26
    IRind{i} = find(IRegion(:)==i);    
    RegPotDNI(i) = sum(AnnPotCellDNI(IRind{i})); %kWh / region / y
end
RegPotDNI(27) =sum(RegPotDNI(1:26));

% Country reg
for i=1:numel(ISOGDP(:,1))
    Cind{i} = find(GISO(:)==ISOGDP(i,1));
    CPotCSP(i) = sum(AnnPotCellDNI(Cind{i})); %kWh / country / y
end
dt = sum(CPotCSP(1:end));
CPotCSP(end+1) =dt;

%% PV Regional and Global Potentials

% Gl_Horiz_NASA{m} % kWh / m2 / day

for r=1:nr
    for c=1:nc
        AnnPotCellPV(r,c) = Gl_Horiz_NASA{13}(r,c) * Area(r,c) * m2tokm2conv * 365; %kWh / cell / y
    end
end

for i=1:26
    RegPotPV(i) = sum(AnnPotCellPV(IRind{i})); % kWh / region / y
end
RegPotPV(27) =sum(RegPotPV(1:26));

for i=1:numel(ISOGDP(:,1))
    CPotPV(i) = sum(AnnPotCellPV(Cind{i})); % kWh / country / y
end
dt = sum(CPotPV(1:end));
CPotPV(end+1) =dt;

%% Insolution kwh/m2/day --> W/m2 conversion

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
    for r=1:nr
        for c=1:nc
            Gl_Horiz_NASA_Insol{m}(r,c) = Gl_Horiz_NASA{m}(r,c) * 41.667;% 1142; % W/m2 insolation
        end
    end
end

%%
% for m=1:13
%     Irr(m) = Gl_Horiz_NASA{m}(101,353);
% end
% plot(Irr)
% imagesc(Gl_Horiz_NASA{13})
