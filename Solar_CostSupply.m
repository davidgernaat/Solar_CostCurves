%% Transform into TIMER curves

for j=1:27
    
    %CSP
    CSP_Pot_cumsum{j} = cumsum(CSP_Pot_IR_S{j}); %kWh
    
    for i=1:numel(CSP_Pot_cumsum{j})
        CSP_Pot_Indexed{j}(i) = CSP_Pot_cumsum{j}(i)/CSP_Pot_cumsum{j}(end); 
    end
    
    %PV
    PV_Pot_cumsum{j} = cumsum(PV_Pot_IR_S{j}); %kWh
    
    for i=1:numel(PV_Pot_cumsum{j})
        PV_Pot_Indexed{j}(i) = PV_Pot_cumsum{j}(i)/PV_Pot_cumsum{j}(end);
    end
    
    %PVres
    PVres_Pot_cumsum{j} = cumsum(PVres_Pot_IR_S{j}); %kWh
    
    for i=1:numel(PVres_Pot_cumsum{j})
        PVres_Pot_Indexed{j}(i) = PVres_Pot_cumsum{j}(i)/PVres_Pot_cumsum{j}(end);
    end
    
    %PVres Urban
    PVres_Pot_cumsumu{j} = cumsum(PVres_Pot_IR_Su{j}); %kWh
    
    for i=1:numel(PVres_Pot_cumsumu{j})
        PVres_Pot_Indexedu{j}(i) = PVres_Pot_cumsumu{j}(i)/PVres_Pot_cumsumu{j}(end);
    end
    
    %PVres Rural
    PVres_Pot_cumsumr{j} = cumsum(PVres_Pot_IR_Sr{j}); %kWh
    
    for i=1:numel(PVres_Pot_cumsumr{j})
        PVres_Pot_Indexedr{j}(i) = PVres_Pot_cumsumr{j}(i)/PVres_Pot_cumsumr{j}(end);
    end
    
    
end

%For PVres paper with 10 regions
for j=1:11
    
    %PVres
    PVres_Pot_cumsum_p{j} = cumsum(PVres_Pot_IR_S_p{j}); %kWh
    
    for i=1:numel(PVres_Pot_cumsum_p{j})
        PVres_Pot_Indexed_p{j}(i) = PVres_Pot_cumsum_p{j}(i)/PVres_Pot_cumsum_p{j}(end);
    end

    
end

%% Find min idx
val = linspace(0,1,100);

for j=1:27
    for i=1:numel(val)
        
        %CSP
        clear tmp
        if isempty(CSP_Pot_Indexed{j})==1; continue; end
        tmp = abs((CSP_Pot_Indexed{j}./val(i))-1);
        [~, idxminCSP{j}(i)] = min(tmp); %index of closest value
        
        %PV
        clear tmp
        if isempty(PV_Pot_Indexed{j})==1; continue; end
        tmp = abs((PV_Pot_Indexed{j}./val(i))-1);
        [~, idxminPV{j}(i)] = min(tmp); %index of closest value
        
        %PVres total
        clear tmp
        if isempty(PVres_Pot_Indexed{j})==1; continue; end
        tmp = abs((PVres_Pot_Indexed{j}./val(i))-1);
        [~, idxminPVres{j}(i)] = min(tmp); %index of closest value
        
        %PVres urban
        clear tmp
        if isempty(PVres_Pot_Indexedu{j})==1; continue; end
        tmp = abs((PVres_Pot_Indexedu{j}./val(i))-1);
        [~, idxminPVresu{j}(i)] = min(tmp); %index of closest value
        
        %PVres rural
        clear tmp
        if isempty(PVres_Pot_Indexedr{j})==1; continue; end
        tmp = abs((PVres_Pot_Indexedr{j}./val(i))-1);
        [~, idxminPVresr{j}(i)] = min(tmp); %index of closest value
    end
end

%For PVres paper with 10 regions
for j=1:11
    for i=1:numel(val)
        %PVres total
        clear tmp
        if isempty(PVres_Pot_Indexed_p{j})==1; continue; end
        tmp = abs((PVres_Pot_Indexed_p{j}./val(i))-1);
        [~, idxminPVres_p{j}(i)] = min(tmp); %index of closest value
    end
end


%% Put in TIMER matrix

a = [0	0.01	0.02	0.03	0.04	0.05	0.06	0.07	0.08	0.09	0.1	0.11	0.12	0.13	0.14	0.15	0.16	0.17	0.18	0.19	0.2	0.21	0.22	0.23	0.24	0.25	0.26	0.27	0.28	0.29	0.3	0.31	0.32	0.33	0.34	0.35	0.36	0.37	0.38	0.39	0.4	0.41	0.42	0.43	0.44	0.45	0.46	0.47	0.48	0.49	0.5	0.51	0.52	0.53	0.54	0.55	0.56	0.57	0.58	0.59	0.6	0.61	0.62	0.63	0.64	0.65	0.66	0.67	0.68	0.69	0.7	0.71	0.72	0.73	0.74	0.75	0.76	0.77	0.78	0.79	0.8	0.81	0.82	0.83	0.84	0.85	0.86	0.87	0.88	0.89	0.9	0.91	0.92	0.93	0.94	0.95	0.96	0.97	0.98	0.99	1];

CostCurveSmthCSP = zeros(100,28);
CostCurveSmthPV = zeros(100,28);
CostCurveSmthPVres = zeros(100,28);
CostCurveSmthPVresu = zeros(100,28);
CostCurveSmthPVresr = zeros(100,28);
CostCurveSmthPVres_p = zeros(100,11);

LFCurveSmthCSP = zeros(100,28);
LFCurveSmthPV = zeros(100,28);
LFCurveSmthPVres = zeros(100,28);
LFCurveSmthPVresu = zeros(100,28);
LFCurveSmthPVresr = zeros(100,28);
LFCurveSmthPVres_p = zeros(100,11);

for i=1:numel(a)-1
    CostCurveSmthCSP(i,1) = a(i);
    CostCurveSmthPV(i,1) = a(i);
    CostCurveSmthPVres(i,1) = a(i);
    CostCurveSmthPVresu(i,1) = a(i);
    CostCurveSmthPVresr(i,1) = a(i);
    CostCurveSmthPVres_p(i,1) = a(i);
    
    LFCurveSmthCSP(i,1) = a(i);
    LFCurveSmthPV(i,1) = a(i);
    LFCurveSmthPVres(i,1) = a(i);
    LFCurveSmthPVresu(i,1) = a(i);
    LFCurveSmthPVresr(i,1) = a(i);
    LFCurveSmthPVres_p(i,1) = a(i);
end

for i=1:27
    COEISort_Indexed_CSP{i} = CSP_COE_IR_S{i}(idxminCSP{i});
    COEISort_Indexed_PV{i} = PV_COE_IR_S{i}(idxminPV{i});
    COEISort_Indexed_PVres{i} = PVres_COE_IR_S{i}(idxminPVres{i});
    COEISort_Indexed_PVresu{i} = PVres_COE_IR_Su{i}(idxminPVresu{i});
    COEISort_Indexed_PVresr{i} = PVres_COE_IR_Sr{i}(idxminPVresr{i});
    
    LFISort_Indexed_CSP{i} = CSP_LF_IR_S{i}(idxminCSP{i});
    LFISort_Indexed_PV{i} = PV_LF_IR_S{i}(idxminPV{i});
    LFISort_Indexed_PVres{i} = PVres_LF_IR_S{i}(idxminPVres{i});
    LFISort_Indexed_PVresu{i} = PVres_LF_IR_Su{i}(idxminPVresu{i});
    LFISort_Indexed_PVresr{i} = PVres_LF_IR_Sr{i}(idxminPVresr{i});
end

%For PVres paper with 10 regions
for i=1:11
    COEISort_Indexed_PVres_p{i} = PVres_COE_IR_S_p{i}(idxminPVres_p{i});
    LFISort_Indexed_PVres_p{i} = PVres_LF_IR_S_p{i}(idxminPVres_p{i});
end

for i=1:27
    % Cost curves
    % CSP
    if isempty(COEISort_Indexed_CSP{i})==1; continue; end
    CostCurveSmthCSP(:,i+1) = COEISort_Indexed_CSP{i};
    
    % PV
    if isempty(COEISort_Indexed_PV{i})==1; continue; end
    CostCurveSmthPV(:,i+1) = COEISort_Indexed_PV{i};
    
    % PVres
    if isempty(COEISort_Indexed_PVres{i})==1; continue; end
    CostCurveSmthPVres(:,i+1) = COEISort_Indexed_PVres{i};
    
    % PVres Urban
    if isempty(COEISort_Indexed_PVresu{i})==1; continue; end
    CostCurveSmthPVresu(:,i+1) = COEISort_Indexed_PVresu{i};
    
    % PVres Rural
    if isempty(COEISort_Indexed_PVresr{i})==1; continue; end
    CostCurveSmthPVresr(:,i+1) = COEISort_Indexed_PVresr{i};
    
    %Load curves 
    % CSP
    if isempty(LFISort_Indexed_CSP{i})==1; continue; end
    LFCurveSmthCSP(:,i+1) = LFISort_Indexed_CSP{i};
    
    % PV
    if isempty(LFISort_Indexed_PV{i})==1; continue; end
    LFCurveSmthPV(:,i+1) = LFISort_Indexed_PV{i};
    
    % PVres
    if isempty(LFISort_Indexed_PVres{i})==1; continue; end
    LFCurveSmthPVres(:,i+1) = LFISort_Indexed_PVres{i};
    
    % PVres urban
    if isempty(LFISort_Indexed_PVresu{i})==1; continue; end
    LFCurveSmthPVresu(:,i+1) = LFISort_Indexed_PVresu{i};
    
    % PVres rural
    if isempty(LFISort_Indexed_PVresr{i})==1; continue; end
    LFCurveSmthPVresr(:,i+1) = LFISort_Indexed_PVresr{i};
end

%For PVres paper with 10 regions
for i=1:11
    % Cost curves
    % PVres
    if isempty(COEISort_Indexed_PVres_p{i})==1; continue; end
    CostCurveSmthPVres_p(:,i+1) = COEISort_Indexed_PVres_p{i};

    %Load curves 
    % PVres
    if isempty(LFISort_Indexed_PVres_p{i})==1; continue; end
    LFCurveSmthPVres_p(:,i+1) = LFISort_Indexed_PVres_p{i};
end

for i=1:27
    % Cost curves
    % CSP
    CostCurveSmthCSP(101,1) = 1;
    CostCurveSmthCSP(101,i+1) = COEISort_Indexed_CSP{i}(end);
    
    % PV
    CostCurveSmthPV(101,1) = 1;
    CostCurveSmthPV(101,i+1) = COEISort_Indexed_PV{i}(end);
    
    % PVres
    CostCurveSmthPVres(101,1) = 1;
    CostCurveSmthPVres(101,i+1) = COEISort_Indexed_PVres{i}(end);
    
    % PVres urban
    CostCurveSmthPVresu(101,1) = 1;
    CostCurveSmthPVresu(101,i+1) = COEISort_Indexed_PVresu{i}(end);
    
    % PVres rural
    CostCurveSmthPVresr(101,1) = 1;
    CostCurveSmthPVresr(101,i+1) = COEISort_Indexed_PVresr{i}(end);
    
    % Load curves 
    % CSP
    LFCurveSmthCSP(101,1) = 1;
    LFCurveSmthCSP(101,i+1) = LFISort_Indexed_CSP{i}(end)-0.05;
    
    % PV
    LFCurveSmthPV(101,1) = 1;
    LFCurveSmthPV(101,i+1) = LFISort_Indexed_PV{i}(end)-0.05;
    
    % PVres
    LFCurveSmthPVres(101,1) = 1;
    LFCurveSmthPVres(101,i+1) = LFISort_Indexed_PVres{i}(end)-0.05;
    
    % PVres urban
    LFCurveSmthPVresu(101,1) = 1;
    LFCurveSmthPVresu(101,i+1) = LFISort_Indexed_PVresu{i}(end)-0.05;
    
    % PVres rural
    LFCurveSmthPVresr(101,1) = 1;
    LFCurveSmthPVresr(101,i+1) = LFISort_Indexed_PVresr{i}(end)-0.05;
end

%For PVres paper with 10 regions
for i=1:11
    % Cost curves    
    % PVres
    CostCurveSmthPVres_p(101,1) = 1;
    CostCurveSmthPVres_p(101,i+1) = COEISort_Indexed_PVres_p{i}(end);
    % Load curves 
    % PVres
    LFCurveSmthPVres_p(101,1) = 1;
    LFCurveSmthPVres_p(101,i+1) = LFISort_Indexed_PVres_p{i}(end)-0.05;
end

%MaxProd
MaxProdCSP = zeros(1,27);
MaxProdPV = zeros(1,27);
MaxProdPVres = zeros(1,27);
MaxProdPVresu = zeros(1,27);
MaxProdPVresr = zeros(1,27);
MaxProdPVres_p = zeros(1,11);

for i=1:27
    MaxProdCSP(:,i) = CSP_Pot_cumsum{i}(end) * 0.0036; %GJ
    MaxProdPV(:,i) = PV_Pot_cumsum{i}(end) * 0.0036; %GJ
    MaxProdPVres(:,i) = PVres_Pot_cumsum{i}(end) * 0.0036; %GJ
    MaxProdPVresu(:,i) = PVres_Pot_cumsumu{i}(end) * 0.0036; %GJ
    MaxProdPVresr(:,i) = PVres_Pot_cumsumr{i}(end) * 0.0036; %GJ
end

%For PVres paper with 10 regions
for i=1:11
    MaxProdPVres_p(:,i) = PVres_Pot_cumsum_p{i}(end) * 0.0036; %GJ
end
    

%% Show
% Costcurves

% h=figure(1);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(CostCurveSmthCSP(1:end-3,i+1));
%     MaxPotential(i) = max(MaxProdCSP(i))/0.0036;
%     maxpot{i} = sprintf('CSP %s (%0.2f PWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('$2010/kWh','FontSize',8)
%     xlabel('Cumulative energy supply (PWh)','FontSize',8)
% end
% 
% h=figure(2);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(CostCurveSmthPV(1:end-3,i+1));
%     MaxPotential(i) = max(MaxProdPV(i))/0.0036;
%     maxpot{i} = sprintf('PV %s (%0.2f PWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('$2010/kWh','FontSize',8)
%     xlabel('Cumulative energy supply (PWh)','FontSize',8)
% end
% 
% h=figure(3);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(CostCurveSmthPVres(1:end-3,i+1));
%     MaxPotential(i) = max(MaxProdPVres(i))/0.0036;
%     maxpot{i} = sprintf('PVres %s (%0.2f TWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('$2010/kWh','FontSize',8)
%     xlabel('Cumulative energy supply (TWh)','FontSize',8)
% end
% 
% h=figure(4);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(CostCurveSmthPVresu(1:end-3,i+1));
%     MaxPotential(i) = max(MaxProdPVresu(i))/0.0036;
%     maxpot{i} = sprintf('PVres Urban %s (%0.2f TWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('$2010/kWh','FontSize',8)
%     xlabel('Cumulative energy supply (TWh)','FontSize',8)
% end
% 
% h=figure(5);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(CostCurveSmthPVresr(1:end-3,i+1));
%     MaxPotential(i) = max(MaxProdPVresr(i))/0.0036;
%     maxpot{i} = sprintf('PVres Rural %s (%0.2f TWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('$2010/kWh','FontSize',8)
%     xlabel('Cumulative energy supply (TWh)','FontSize',8)
% end

%Loadfactor
% h=figure(1);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(LFCurveSmthCSP(:,i+1));
%     MaxPotential(i) = max(MaxProdCSP(i))/0.0036;
%     maxpot{i} = sprintf('CSP %s (%0.2f PWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('Load factor','FontSize',8)
%     xlabel('Cumulative energy supply (PWh)','FontSize',8)
% end
% 
% h=figure(2);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(LFCurveSmthPV(:,i+1));
%     MaxPotential(i) = max(MaxProdPV(i))/0.0036;
%     maxpot{i} = sprintf('PV %s (%0.2f PWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('Load factor','FontSize',8)
%     xlabel('Cumulative energy supply (PWh)','FontSize',8)
% end
% 
% h=figure(3);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(LFCurveSmthPVres(:,i+1));
%     MaxPotential(i) = max(MaxProdPVres(i))/0.0036;
%     maxpot{i} = sprintf('PVres %s (%0.2f TWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('Load factor','FontSize',8)
%     xlabel('Cumulative energy supply (TWh)','FontSize',8)
% end
% 
% h=figure(4);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(LFCurveSmthPVresu(:,i+1));
%     MaxPotential(i) = max(MaxProdPVresu(i))/0.0036;
%     maxpot{i} = sprintf('PVres urban %s (%0.2f TWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('Load factor','FontSize',8)
%     xlabel('Cumulative energy supply (TWh)','FontSize',8)
% end
% 
% h=figure(5);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(LFCurveSmthPVresr(:,i+1));
%     MaxPotential(i) = max(MaxProdPVresr(i))/0.0036;
%     maxpot{i} = sprintf('PVres rural %s (%0.2f TWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('Load factor','FontSize',8)
%     xlabel('Cumulative energy supply (TWh)','FontSize',8)
% end

%% Save
if setoutput==true;

    pathname = fileparts(root);
    
    % CostCurves
    
    % CSP
    if setCSPPL==0
        file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\CostCurveSmthCSP.dat'));
    elseif setCSPPL==1
        if multifactor==0
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\CostCurveSmthCSP_PL.dat'));
        elseif multifactor==1
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\CostCurveSmthCSP_PLDL.dat'));
        end
    end
    txt=sprintf('real main.CSPCostCurveSmth[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,CostCurveSmthCSP,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % PV
    if setPVPL==0
        file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\CostCurveSmthPV.dat'));
    elseif setPVPL==1
        if multifactor==0
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\CostCurveSmthPV_PL.dat'));
        elseif multifactor==1
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\CostCurveSmthPV_PLDL.dat'));
        end
    end
    txt=sprintf('real main.PVCostCurveSmth[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,CostCurveSmthPV,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % PVres
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\CostCurveSmthPVres.dat'));
    txt=sprintf('real main.PVresCostCurveSmth[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,CostCurveSmthPVres,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % PVres urban
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\CostCurveSmthPVresUrb.dat'));
    txt=sprintf('real main.PVresCostCurveSmthUrb[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,CostCurveSmthPVresu,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % PVres Rural
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\CostCurveSmthPVresRur.dat'));
    txt=sprintf('real main.PVresCostCurveSmthRur[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,CostCurveSmthPVresr,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % PVres for paper
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\CostCurveSmthPVres_p.dat'));
    txt=sprintf('real main.PVresCostCurveSmth[11](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,CostCurveSmthPVres_p,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    
    % MaxProd
    
    %CSP
    if setCSPPL==0
        file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\MaxProdCSP.dat'));
    elseif setCSPPL==1
        if multifactor==0
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\MaxProdCSP_PL.dat'));
        elseif multifactor==1
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\MaxProdCSP_PLDL.dat'));
        end
    end
    txt=sprintf('real main.CSP_TechPotRegion[27] = ');
    dlmwrite(file,txt,'');
    dlmwrite(file,MaxProdCSP,'-append','delimiter','\t');
    dlmwrite(file,';','-append');
    
    % PV
    if setPVPL==0
        file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\MaxProdPV.dat'));
    elseif setPVPL==1
        if multifactor==0
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\MaxProdPV_PL.dat'));
        elseif multifactor==1
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\MaxProdPV_PLDL.dat'));
        end
    end
    txt=sprintf('real main.PV_TechPotRegion[27] = ');
    dlmwrite(file,txt,'');
    dlmwrite(file,MaxProdPV,'-append','delimiter','\t');
    dlmwrite(file,';','-append');
    
    % PVres
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\MaxProdPVres.dat'));
    txt=sprintf('real main.PVres_TechPotRegion[27] = ');
    dlmwrite(file,txt,'');
    dlmwrite(file,MaxProdPVres,'-append','delimiter','\t');
    dlmwrite(file,';','-append');
    
    % PVres urban
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\MaxProdPVresUrb.dat'));
    txt=sprintf('real main.PVres_TechPotRegionUrb[27] = ');
    dlmwrite(file,txt,'');
    dlmwrite(file,MaxProdPVresu,'-append','delimiter','\t');
    dlmwrite(file,';','-append');
    
    % PVres rural
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\MaxProdPVresRur.dat'));
    txt=sprintf('real main.PVres_TechPotRegionRur[27] = ');
    dlmwrite(file,txt,'');
    dlmwrite(file,MaxProdPVresr,'-append','delimiter','\t');
    dlmwrite(file,';','-append');
    
    % PVres for paper
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\MaxProdPVres_p.dat'));
    txt=sprintf('real main.PVres_TechPotRegion[11] = ');
    dlmwrite(file,txt,'');
    dlmwrite(file,MaxProdPVres_p,'-append','delimiter','\t');
    dlmwrite(file,';','-append');
    
    % Load factor decline curves
    
    %CSP
    if setCSPPL==0
        file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\LoadCurveSmthCSP.dat'));
    elseif setCSPPL==1
        if multifactor==0
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\LoadCurveSmthCSP_PL.dat'));
        elseif multifactor==1
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\LoadCurveSmthCSP_PLDL.dat'));
        end
    end
    txt=sprintf('real main.CSPLoadCurveSmth[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,LFCurveSmthCSP,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % PV
    if setPVPL==0
        file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\LoadCurveSmthPV.dat'));
    elseif setPVPL==1
        if multifactor==0
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\LoadCurveSmthPV_PL.dat'));
        elseif multifactor==1
            file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\LoadCurveSmthPV_PLDL.dat'));
        end
    end
    txt=sprintf('real main.PVLoadCurveSmth[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,LFCurveSmthPV,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % PVres
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\LoadCurveSmthPVres.dat')); 
    txt=sprintf('real main.PVresLoadCurveSmth[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,LFCurveSmthPVres,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % PVres urban
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\LoadCurveSmthPVresUrb.dat')); 
    txt=sprintf('real main.PVresLoadCurveSmthUrb[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,LFCurveSmthPVresu,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % PVres rural
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\LoadCurveSmthPVresRur.dat')); 
    txt=sprintf('real main.PVresLoadCurveSmthRur[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,LFCurveSmthPVresr,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % PVres for paper with 10 regions
    file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\LoadCurveSmthPVres_p.dat')); 
    txt=sprintf('real main.PVresLoadCurveSmth[11](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,LFCurveSmthPVres_p,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
end