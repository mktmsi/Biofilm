clear
CheckV=1;
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/Syuron/Single/1/Results/Good/ArrayZ';
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/Syuron/InPhase/Results/Good/ArrayV';
folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/AfterSyuron/Inter/paratan/Results/ParameterSet15/ArrayU';
folderNameV='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/AfterSyuron/Inter/paratan/Results/ParameterSet15/ArrayV';
folderNameSignal='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/AfterSyuron/Inter/paratan/Results/ParameterSet15/SignalFlags';

cd(folderName);
%dir('*.csv')
MyFolderInfo = dir('*.csv');
NumFiles=length(MyFolderInfo)


for i=1:NumFiles-1

if CheckV==1
    
    fileNameSignaling=strcat(folderNameSignal,'/',num2str(i),'.csv');
    fileSignaling=csvread(fileNameSignaling);
    [row,col]=size(fileSignaling);
    fileSignaling=100*fileSignaling;
    max(max(fileSignaling))
    fileSignaling= ones('like',fileSignaling)+fileSignaling;
    
    %clf
    fileNameV=strcat(folderNameV,'/',num2str(i),'.csv');
    fileDataV=csvread(fileNameV);
    [row,col]=size(fileDataV);
    [XX,YY]=meshgrid(1:col,1:row);
    fileDataV(fileDataV<0.0)=0.0;
    fileDataV(fileDataV>0.01)=50.0;
    fileDataV=fileDataV.*fileSignaling;
    %fileDataV(fileDataV>0.0 & fileDataV<=22.0)=100.0;
end

%cd(folderName);
%dir('*.csv')
%MyFolderInfo = dir('*.csv');
%NumFiles=length(MyFolderInfo)

    fileName=strcat(folderName,'/',num2str(i),'.csv');
    i;
    fileData=csvread(fileName);
    [row,col]=size(fileData);
    [XX,YY]=meshgrid(1:col,1:row);
    size(fileData);
    size(XX);
    s=surf(XX,YY,fileData,'FaceAlpha',0.5);
    s.CData=fileDataV;
    colorbar
    %surf(XX,YY,fileData)
    zlim([0 500])
    xlim([0 col-1]) 
    ylim([0 col-1]) 
    xlabel('x')
    ylabel('y')
    zlabel('v')
    %pause(0.01);
    F(i+1) = getframe(gcf);
    
end

% write to video
v = VideoWriter('V.avi');
v.FrameRate = 30; % Framerate

open(v);
writeVideo(v,F);
close(v);