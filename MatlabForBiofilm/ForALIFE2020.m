clear
FlagSingle=0;
rowEnd=1000;
%%グラフの制御
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=5;
%XlimStart=2.0;
%XlimEnd=3.1;
folderName='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/1D/Working/Single/1/Results/ParameterSet0';
cd(folderName);
fileName=strcat(folderName,'/SumV.csv');
%fileName=strcat(folderName,'/V_para.csv');
%fileName=strcat(folderName,'/SumV.csv');
fileData=csvread(fileName);
[row,col]=size(fileData);
if row==282
    fileData=fileData(1:row-1,:);
end
[row,col]=size(fileData);
%%%2つのバイオフィルムのサイズの時間推移を表示%%%%%%%%%%%
if FlagSingle==0
    start=1;
    dend=row;
    figure;
    b=cast((col-2)/2,'int64');
    %DataLeft=sum(fileData(start:dend,2:b),2);
    % DataRight=sum(fileData(start:dend,b:col),2);
    % DataLeft=DataLeft.';
    % DataRight=DataRight.';
    DataLeft=fileData(start:dend,5);
    %DataRight=fileData(start:dend,3);
     XData=(fileData(start:dend,1)/1000000);%ALIFE2020
   % XData=(fileData(start:dend,1)/1000);%ALIFE2020
    plot(XData,DataLeft,'LineWidth',LineHaba);
    hold on;
    %plot(XData,DataRight,'LineWidth',LineHaba);
    %c=legend('Biofilm1','Biofilm2');%凡例を制御
    c.FontSize=HanreiSize;
    c.Location='northwest';
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%軸の太さ
        XlimStart=4.90;%for alife
   % XlimEnd=5.00;%for alife
   % xlim([XlimStart XlimEnd]);% for alife
    % xlim([min(XData) inf]);
    box off%余分な右のy軸と上のx軸を削除
    %軸ラベルの追加
    ax=gca;
    xlabel('Time [\times 10^6 steps]','FontSize',LabelFontSize,'FontWeight','bold')
    ylabel('Biofilm size','FontSize',LabelFontSize,'FontWeight','bold')
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba)%メモリのない上側のx軸と右側のy軸を追加

    %グラフ保存

    filenameBlend=strcat(folderName,'/Expansion.png');
    saveas(gcf,filenameBlend);
    %グラフの微分値を取得
    dDataLeft=diff(DataLeft);
    %dDataRight=diff(DataRight);
    [row2,col2]=size(dDataLeft);  
    XDataFiltered=XData(1:row2)
    figure;
    %%左側のバイオフィルム
    ax=plot(XDataFiltered,dDataLeft,'LineWidth',LineHaba);
    %findpeaks(FilterdDataLeft,XDataFiltered)
    hold on;
    %%右側のバイオフィルム をプロット
    %plot(XDataFiltered,dDataRight,'-.','LineWidth',LineHaba);
    c=legend('Biofilm1','Biofilm2');%凡例を制御
    c.FontSize=HanreiSize;
    c.Location='northwest';
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%軸の太さ

    %xlim([-5 XlimEnd]);
    % xlim([min(XData) inf]);
    box off%余分な右のy軸と上のx軸を削除
    %軸ラベルの追加
    ax=gca;
    XlimStart=1.5;
    XlimEnd=3.0;
    xlim([XlimStart XlimEnd]);% for alife
    %xlabel('Time [\times 10^3 steps]','FontSize',LabelFontSize,'FontWeight','bold')
    ylabel('Expansion rate ','FontSize',LabelFontSize,'FontWeight','bold')
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba)%メモリのない上側のx軸と右側のy軸を追加
    
    filenameBlend=strcat(folderName,'/dExpansionNotFiltered.png');
    saveas(ax,filenameBlend);
end