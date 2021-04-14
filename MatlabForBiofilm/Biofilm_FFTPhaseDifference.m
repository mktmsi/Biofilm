clear;

FlagSingle=0;
rowEnd=1000;
%%�O���t�̐���
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=5;
XlimStart=2.0;
XlimEnd=3.1;


folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/FullPaper/AutoModes/Results0609/ParameterSet';
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/ALIFE2020CameraReady/AntiPhase';
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/OldForALIFE2020/SmallField/InPhase';
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/ParameterSet3';
folderName=strcat(folderName,num2str(21));
cd(folderName);
fileName=strcat(folderName,'/Vsum.csv');
%fileName=strcat(folderName,'/V_para.csv');
%fileName=strcat(folderName,'/SumV.csv');
fileData=csvread(fileName);
[row,col]=size(fileData);
%%%2�̃o�C�I�t�B�����̃T�C�Y�̎��Ԑ��ڂ�\��%%%%%%%%%%%
if FlagSingle==0
    start=2;
    dend=row;
    figure;
    b=cast((col-2)/2,'int64');
    %DataLeft=sum(fileData(start:dend,2:b),2);
    % DataRight=sum(fileData(start:dend,b:col),2);
    % DataLeft=DataLeft.';
    % DataRight=DataRight.';
    DataLeft=fileData(start:dend,2);
    DataRight=fileData(start:dend,3);
    %XData=(fileData(start:dend,1)/1000000);%ALIFE2020
    XData=(fileData(start:dend,1)/10000);%Fullpaper
    plot(XData,DataLeft,'LineWidth',LineHaba);
    hold on;
    plot(XData,DataRight,'LineWidth',LineHaba);
    c=legend('Biofilm1','Biofilm2');%�}��𐧌�
    c.FontSize=HanreiSize;
    c.Location='northwest';
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%���̑���
    XlimStart=max(XData)/2;
    XlimEnd=max(XData);
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
    dDataLeft=diff(DataLeft);
    dDataRight=diff(DataRight);
    [row2,col2]=size(dDataLeft);
    %plot(dDataLeft)
    
    %%%%%%���[�p�X�t�B���^�ɒʂ����g�嗦��\��%%%%%%%
    %���[�p�X�t�B���^ �̐݌v
    fc =0.40;
    Wn = (2.0/row2)*fc;%�J�b�g�I�t
    % t = linspace(0,1,row2);
    %b = fir1(200,Wn,'low');%2D
    b = fir1(30,Wn,'low');%1D
    %�t�B���^ ��K�p
    FilterdDataLeft=filter(b,1,dDataLeft);
    FilterdDataRight=filter(b,1,dDataRight);
    [rowFilterd,colFiltered]=size(FilterdDataLeft);
    XDataFiltered=XData(1:rowFilterd);
    
    fs = 100;
    t = 0:1/fs:1-1/fs;
    XlimStart=1;
    XlimEnd=280;
    %x=cos(2*pi*15*t - pi/4) + cos(2*pi*15*t);
    %figure;
    %plot(t,x);
    x = DataRight(XlimStart:XlimEnd);
    %x = DataLeft(XlimStart:XlimEnd)-DataRight(XlimStart:XlimEnd);
    figure;
    y = fft(x);
    z = fftshift(y);
   
    ly = length(y);
    f = (-ly/2:ly/2-1)/ly*fs;
  
    figure;
    stem(f,abs(z))
    xlabel 'Frequency (Hz)'
    ylabel '|y|'
    grid
    
    figure;
    tol = 1e+8;
    z(abs(z) < tol) = 0;
    
    theta = angle(z)
    
    stem(f,theta/pi)
    xlabel 'Frequency (Hz)'
    ylabel 'Phase / \pi'
    grid
    
    
    
    
    
    %%%%�o�C�I�t�B�����Ԃ̈ʑ��������߂�
    
end



