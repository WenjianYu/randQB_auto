%% read image
imname= 'image1.jpg';
X= imread(imname);
X= double(X)/255;
if ndims(X) == 3
   image(X)
   [m,n,p]= size(X);
   if m>n,
       A= [X(:,:,1), X(:,:,2), X(:,:,3)]; 
   else
       A= [X(:,:,1); X(:,:,2); X(:,:,3)]; 
   end
else
   image(255*X);
   colormap(gray(256));
   A= X;
end
axis image
axis off
drawnow
clear X;



