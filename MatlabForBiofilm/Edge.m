%filenameV = '/Users/mikamitaishi/Desktop/SpatioTemporalPlot/ver2/v_values.csv';
filenameV = '/Users/mikamitaishi/Desktop/SpatioTemporalPlotOriginal/v_values.csv';
filenameU = '/Users/mikamitaishi/Desktop/SpatioTemporalPlotOriginal/u_values.csv';
V = csvread(filenameV);
U = csvread(filenameU);
%V = transpose(V);

V = V(1:600,1:900);
U = U(1:600,1:900);
size(V)
%V1 = V(:);
%V=V1;
%N(:,:,1)=M
normV = V ./ max(V(:));
%imU=imagesc(normV);
normU = U ./ max(U(:));
thrash=0.125;
normV=normV>thrash;
Vmax=thrash*max(V(:));
%N=zeros(size1,size2,3)
%N(:,:,2)=normM

%imV=imagesc(normV);
%hold on

%エッジ検出
[xsize,ysize]=size(normV);
filter=(1/8)*[1 2 1; 0 0 0; -1 -2 -1];%sobel 縦
for i=2:xsize-1
    for j=2:ysize-1 
        temp=normV(i-1:i+1, j-1:j+1); 
        result=temp.*filter; 
        ImgOut(i, j)=sum(result(:));
    end
end
ImgOut=abs(ImgOut);
%ImgOut = ImgOut(1:599,1:900);
%normU=normU(1:599,1:900);
%imshow(ImgOut)
size(ImgOut)
size(normU)
%エッジ座標の検出
[row,col] = find(ImgOut);
sizeRow=size(row)

[sizex,sizey]=size(normU)
X=zeros(sizex,sizey);
for i=1:sizeRow
        X(row(i,1),col(i,1))=1.0;
end
imV=imshow(X)
saveas(imV,'/Users/mikamitaishi/Desktop/SpatioTemporalPlot/SpatioTemporalV.png')
saveas(imU,'/Users/mikamitaishi/Desktop/SpatioTemporalPlot/SpatioTemporalU.png')
imV=imread('/Users/mikamitaishi/Desktop/SpatioTemporalPlot/SpatioTemporalV.png')
%X=zeros(,)
colormap('parula')
imU=imagesc(normU);
hold on

%
colorbar 





