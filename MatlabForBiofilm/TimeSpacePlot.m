clear
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/AfterSyuron/Inter/paratan/Results/B-3';
folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/Syuron/NewResults/Single';
filenameV = strcat(folderName,'/V_para9.csv')
filenameU = strcat(folderName,'/U_para9.csv');
filenameZ = strcat(folderName,'/Z_para9.csv');
V = csvread(filenameV);
size(V)
U = csvread(filenameU);
size(U);
Z = csvread(filenameZ);

LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=3;


%データの整形
Sizes=size(Z);
gyou_end=1000;
retu_first=2;
retu_end=102;
TimeAxis=V(1:gyou_end,1);
V = V(1:gyou_end,retu_first:retu_end);
thre=1.0;
V(V>thre)=thre;%???
%V(V>0.1)=1;
%V(V~=1)=0;
%V=transpose(V);
%エッジ検出
BW2 = edge(V);

U = U(1:gyou_end,retu_first:retu_end);
Umax=max(U(:))
normU = U / max(U(:));
%normU=1 - normU;
%normU=transpose(normU);
%normU(normU<0.8)=0;

Z = Z(1:gyou_end,retu_first:retu_end);
Z(Z>10.0)=10.0;
%normZ = Z / max(Z(:));
%normZ=1 - normZ;
%%%%%%%%%%画像表示%%%%%%%%%%%%
clf
figure
figure('Name','V_value')
%subplot(2,1,1)
imV=imagesc(V);
colormap parula
colorbar

%clf
%figure
%figure('Name','Z_value')
%subplot(2,1,1)
%imV=imagesc(Z)
%colormap parula
%colorbar

%figure('Name','V_edges')
%hold on
%imV2=imagesc(BW2)
BW2=1-BW2;
%colormap gray
%saveas(imV2,'/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/BiofilmAuto/1/Results/ParameterSet0/VEdge.png')

%{
figure
subplot(2,1,1)
%figure('Name','U_values')
imU=imagesc(normU)
colormap parula
%}
%saveas(im,'/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/BiofilmAuto/1/Results/ParameterSet0/Uvalue.png')

%%%%%Vのエッジとエサ量との重ね合わせ画像
figure('Name','Blend')
%subplot(2,1,2)
Blend=immultiply(BW2,normU);
Blend = (rot90(Blend));
imBl=imagesc(Blend)
filenameBlend=strcat(folderName,'/Blend1.png');
TimeAxis=TimeAxis/1000;
%yticklabels(TimeAxis);
colormap parula
cmin=0.3;
cmax=0.8;
caxis([cmin cmax])%色の範囲
c = colorbar;
c.Ticks=([cmin,(cmin+cmax)/2,cmax]);
c.TickLabels=({'0','30','60'});
minT=min(TimeAxis);
maxT=max(TimeAxis)
xticks([1 :220 :gyou_end+1]);
xticklabels({'0','10','20','30','40'});
set(gca,'FontSize',MemoriSize);
yticks([1,50,100])
yticklabels({'-r_0','0','r_0'})
%xticks({'-r_0','0','r_0'})
xlabel('Time [\times 10^3 steps]','FontSize',LabelFontSize,'FontWeight','bold')
ylabel('Coordinate x','FontSize',LabelFontSize,'FontWeight','bold')
ax=gca;
box off%余分な右のy軸と上のx軸を削除
set(gca,'LineWidth',LineHaba);
axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba);%目盛りのない上側のx軸と右側のy軸を追加

saveas(imBl,filenameBlend);


