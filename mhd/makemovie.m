function [ A ] = makemovie( data, numframes, clims)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

fig1 = figure(1)

winsize = get(fig1,'Position');
winsize(1:2) = [0,0];

A=moviein(numframes,fig1,winsize);
for i=1:numframes
   image(:,:) = (data(i,:,:));
   imagesc(image,clims)
   axis xy
   colorbar
  
   A(:,i)=getframe(fig1,winsize);
 end

end

