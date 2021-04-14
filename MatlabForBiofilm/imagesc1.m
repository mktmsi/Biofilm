% script for Speed and CoT
% author:   Kyoichi Akiyama
% date:     25-01-2019
% Nan,inf->zero
close all;
clear all;

%% load file
filename = 'control1_6.csv';
A = csvread(filename,1,0);
bins = 10;
size = 20;

%frequency
freq0 = A(1:bins,5);%A(1:last,5th)
freq = freq0.';

%amplitude
amp0 = A(:,4);
amp1 = reshape(amp0,[bins,bins]);
amp = fliplr(amp1(1,:));%opposite order
% amp = amp1(1,:);%opposite order

%speed
speed0 = A(:,2);
speed1 = reshape(speed0,[bins,bins]);
speed  = flipud(speed1.');

%CoT
cot0 = A(:,3);
cot1 = reshape(cot0,[bins,bins]);
cot2  = flipud(cot1.');
cot = 1./cot2;
% cot2(cot2>10) = 0;%minus ->0
%% each graph
figure;
speedlims = [0 3];
imagesc(freq, amp(1:8), speed(1:8,:), speedlims);
% imagesc(freq, amp, speed);
colorbar;
colormap('hot');
set(gca,'YDir','normal');
% title('Speed [m/sec]');
xlabel('\omega');
ylabel('{\it C_{amp}}');
set(gca,'YTick',[0.04:0.08:0.2]);
set(gca,'XTick',[0.5 1.0 1.5]);
set(gca,'xticklabels',{'\pi','2\pi','3\pi'});
set(gca,'FontSize',size);
set(gca,'FontName','times new roman'); 

figure;
cotlims = [0.2 0.7];
imagesc(freq, amp(1:8), cot(1:8,:), cotlims);
% imagesc(freq, amp, cot2);
colorbar;
colormap('cool');
set(gca,'YDir','normal');
% title('CoT^{-1}');
xlabel('\omega');
ylabel('{\it C_{amp}}');
set(gca,'YTick',[0.04:0.08:0.2]);
set(gca,'XTick',[0.5 1.0 1.5]);
set(gca,'xticklabels',{'\pi','2\pi','3\pi'});
set(gca,'FontSize',size);
set(gca,'FontName','times new roman'); 
