clear

FlagSingle=0;
rowEnd=1000;
%%�O���t�̐���
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=5;
%XlimStart=2.0;
%XlimEnd=3.1;

folderName='/Users/mikamihiroshi/Desktop/Single/Results/ForPoster';
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/ALIFE2020CameraReady/AntiPhase';
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/OldForALIFE2020/SmallField/InPhase';
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/ParameterSet3';
%folderName=strcat(folderName,num2str(79));
%cd(folderName);
%fileName=strcat(folderName,'/Vsum.csv');
%fileName=strcat(folderName,'/V_para.csv');
fileName=strcat(folderName,'/SumV.csv');
fileData=csvread(fileName);
[row,col]=size(fileData);
%%%2�̃o�C�I�t�B�����̃T�C�Y�̎��Ԑ��ڂ�\��%%%%%%%%%%%
if FlagSingle==0
    start=1;
    dend=row;
    figure;
    b=cast((col-2)/2,'int64');
    DataLeft=fileData(start:dend,2);
    XData=(fileData(start:dend,1)/1000000);%ALIFE2020
    %XData=(fileData(start:dend,1)/10000);%Fullpaper
    plot(XData,DataLeft,'LineWidth',LineHaba);
    hold on;
    c=legend('Biofilm1','Biofilm2');%�}��𐧌�
    c.FontSize=HanreiSize;
    c.Location='northwest';
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%���̑���
    XlimStart=4.0;
    XlimEnd=4.20;
    xlim([XlimStart XlimEnd]);
    % xlim([min(XData) inf]);
    box off%�]���ȉE��y���Ə��x�����폜
    %�����x���̒ǉ�
    ax=gca;
    %ax.XTick=[0 200 400 600 800 1000];%�����I�Ɏ��̍��݂�200���ɂ���
    %axis([0 1000 0 inf])
    %xticklabels({'0','20','40','60','80','100'})
    xlabel('Time [\times 10^3 steps]','FontSize',LabelFontSize,'FontWeight','bold')
    ylabel('Biofilm size','FontSize',LabelFontSize,'FontWeight','bold')
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba)%�������̂Ȃ��㑤��x���ƉE����y����ǉ�
    %�O���t�ۑ�
    filenameBlend=strcat(folderName,'/Expansion.png');
    saveas(gcf,filenameBlend);
    %�O���t�̔����l���擾
    dDataLeft={};
    dDataLeft=diff(DataLeft);
    [row2,col2]=size(dDataLeft);
    figure;
    %%�����̃o�C�I�t�B����
    XData=XData(1:row2)
    ax=plot(XData,dDataLeft,'LineWidth',LineHaba);
    %findpeaks(FilterdDataLeft,XDataFiltered)
    c=legend('Biofilm1');%�}��𐧌�
    c.FontSize=HanreiSize;
    c.Location='northwest';
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%���̑���

    xlim([XlimStart XlimEnd]);
    % xlim([min(XData) inf]);
    box off%�]���ȉE��y���Ə��x�����폜
    %�����x���̒ǉ�
    ax=gca;
    xlabel('Time [\times 10^6 steps]','FontSize',LabelFontSize,'FontWeight','bold')
    ylabel('Expansion rate ','FontSize',LabelFontSize,'FontWeight','bold')
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba)%�������̂Ȃ��㑤��x���ƉE����y����ǉ�
   
    filenameBlend=strcat(folderName,'/dExpansion.png');
    %%�O���t��ۑ�
    saveas(ax,filenameBlend);

end


%%%�P�o�C�I�t�B����%%%%
if FlagSingle==1
    Data=fileData(1:rowEnd,2);
    XData=(fileData(1:rowEnd,1)/1000);
    f=figure;
    plot(XData,Data,'LineWidth',LineHaba);
    xlim([min(XData) max(XData)]);
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%���̑���
    box off%�]���ȉE��y���Ə��x�����폜
    %�����x���̐���
    ax=gca;
    %ax.XTick=[0 100 200 300 400 500 600 700]%�����I�Ɏ��̍��݂�100���ɂ���
    %axis([0 700 0 500])
    %xticklabels({'0','10','20','30','40','50','60','70'})%�ڐ���̖��O�t������
    %xticklabels(strXData);
    xlabel('Time [\times 10^3 steps]','FontSize',LabelFontSize,'FontWeight','bold');
    ylabel('Biofilm size','FontSize',LabelFontSize,'FontWeight','bold');
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba);%�ڐ���̂Ȃ��㑤��x���ƉE����y����ǉ�
    %�O���t�ۑ�
    filenameBlend=strcat(folderName,'/Expansion.png');
    saveas(gcf,filenameBlend);
    %�O���t�̔����l���擾
    dData=diff(Data);
    [row2,col2]=size(dData);
    %plot(dDataLeft)
    %%%%%%���[�p�X�t�B���^�ɒʂ����g�嗦��\��%%%%%%%
    %���[�p�X�t�B���^ �̐݌v
    fc =10000.0;
    Wn = (2.0/row2)*fc;
    t = linspace(0,0.0001,row2);
    b = fir1(200,Wn,'low');
    %�t�B���^ ��K�p
    %dData=dData(450:925);
    FilterdDataLeft=filter(b,1,dData);
    size(FilterdDataLeft);
    %�O���t�\��
    [row3,col3]=size(dData);
    %XData=XData(450:925);
    max(XData);
    figure;
    ax=plot(XData,FilterdDataLeft,'LineWidth',LineHaba);
    xlim([min(XData) max(XData)]);
    %xticklabels({'0','20','40','60','80','100'})
    %xticklabels({'70','75','80','85','90','95','100'})
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);
    xlabel('Time [\times 10^3 steps]','FontSize',LabelFontSize,'FontWeight','bold');
    ylabel('Expanding rate [/step]','FontSize',LabelFontSize,'FontWeight','bold');
    ax=gca;
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba);
    filenameBlend=strcat(folderName,'/dExpansion.png');
    saveas(gcf,filenameBlend);
end

