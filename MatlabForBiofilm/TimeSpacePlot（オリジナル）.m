%filenameV = '/Users/mikamitaishi/Desktop/SpatioTemporalPlot/v_values.csv';
filenameV = '/Users/mikamitaishi/Desktop/forALIFE2019_iMac6/v_values.csv';
V = csvread(filenameV);
%N(:,:,1)=M
normV = V ./ max(V(:));
%N=zeros(size1,size2,3)
%N(:,:,2)=normM

colormap('parula')
imV=imagesc(normV)
caxis([0.0 2.0])
colorbar 
%saveas(imV,'/Users/mikamitaishi/Desktop/SpatioTemporalPlot/SpatioYemporalV_12.png')



