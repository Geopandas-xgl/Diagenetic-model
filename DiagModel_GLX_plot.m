function DiagModel_GLX_plot
%% ============================ Load data =================================
clear all
warning off
load("data\ModelData.mat");
data = readtable("data\Tingri_dC_dO.xlsx","Sheet","C and O Isotope");
%% ==================== d13C phases for Tingri section =====================

PhaseDepth = [4,12,16.2,16.9,17.5,20.8,21.6];
[Phaserow,~,~] = find(data.Depth == PhaseDepth);
ph = array2table([Phaserow';PhaseDepth],"VariableNames",["a","b","c","d","e","f","g"]);
%% ==================== Select specified time and box =====================
box_sel = [1,5,10,20,50];
t_sel = [0.05e6,0.1e6,0.15e6,0.25e6,0.5e6];
t_idx_sel = [2,length(ModelData(1).WR)];

for N = 1:2 
    for i = 1:length(t_sel)
        [~,idx]= min(abs(ModelData(N).WR(:,1)-t_sel(i))); 
        t_idx_sel(N,i) = idx; 
    end
end
%% ========================= Figure =======================================
figure("Position",[0,0,350,300])
index_a = find(data.Depth ==ph.a(2));
index_c = find(data.Depth ==ph.c(2));
for j = t_idx_sel(1,:)
    plot(ModelData(1).Solid.dO(j,:), ModelData(1).Solid.dC(j,:),"LineStyle","--","Color",[0 0 0],'LineWidth',0.8); hold on   
end
plot(ModelData(1).Solid.dO(j,:), ModelData(1).Solid.dC(j,:),"Color",[0 0.4470 0.7410],'LineWidth',1.2);
for j = 1:length(box_sel)
    plot(ModelData(1).Solid.dO(:,box_sel(j)),ModelData(1).Solid.dC(:,box_sel(j)),"Color",[0 0 0],'LineWidth',0.8); hold on
end

for j = t_idx_sel(2,:)
    plot(ModelData(2).Solid.dO(j,:), ModelData(2).Solid.dC(j,:),"LineStyle","--","Color",[0.7 0.7 0.7],'LineWidth',0.8); hold on   
end
plot(ModelData(2).Solid.dO(j,:), ModelData(2).Solid.dC(j,:),"Color",[0.7 0.7 0.7],'LineWidth',1.2);
for j = 1:length(box_sel)
    plot(ModelData(2).Solid.dO(:,box_sel(j)),ModelData(2).Solid.dC(:,box_sel(j)),"Color",[0.7 0.7 0.7],'LineWidth',0.8); hold on
end
sct1 = scatter(data.dO(1:index_a-1),data.dC(1:index_a-1), 20, "Marker","o","MarkerEdgeColor","#000000", "MarkerFaceColor","#FF9D3B","LineWidth",0.5); hold on
sct2 = scatter(data.dO(index_a:index_c-1),data.dC(index_a:index_c-1), 20, "Marker","o","MarkerEdgeColor","#000000", "MarkerFaceColor","#EC1C24","LineWidth",0.5); 
sct3 = scatter(data.dO(index_c:end),data.dC(index_c:end), 20, "Marker","o","MarkerEdgeColor","#000000", "MarkerFaceColor","#00ADEE","LineWidth",0.5); 
sct4 = scatter(data.dO(end-1:end),data.dC(end-1:end), 20, "Marker","o","MarkerEdgeColor","#000000", "MarkerFaceColor","#DFDFE1","LineWidth",0.5);
uistack(sct1,"up");
xlim([-18,1]);
xticks(-15:5:0);
ticklabels("x",xticks,1,"F");
xlabel("\delta^{18}O (‰)");
ylim([-16,5]);
yticks(-15:5:5);
ylabel("\delta^{13}C (‰)");
[lgd,icons]= legend([sct1,sct2,sct3,sct4],["Bottom","Pre-CIE","Main CIE","Top"],"Location","southeast","FontSize",6,"LineWidth",0.5,"EdgeColor","#EEEEEF","FontName","Times New Roman");
lgd.ItemTokenSize = [20,20];
icons = findobj(icons,'Type','patch');
icons = findobj(icons,'Marker','none','-xor');
set(icons,'MarkerSize',5); 
set(gca,"linewidth", 1,"FontSize",12,"FontName", "Times New Roman","TickLength",[0.02,0.02],"Layer","top");
end