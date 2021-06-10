% simple code to test registration by |fhat| in 2d ;

G4 = imread('test_G4.png'); G5 = imread('test_G5.png');
Ax = (255-mean(G4,3))/255; Bx = (255-mean(G5,3))/255;
Ak = recenter2(fft2(recenter2(Ax)));
Bk = recenter2(fft2(recenter2(Bx)));
disp_flag=0;
if (disp_flag);
subplot(2,3,1+0*3); imagesc(Ax);
subplot(2,3,2+0*3); imagesc(abs(Ak));
subplot(2,3,3+0*3); imagesc(angle(Ak));
subplot(2,3,1+1*3); imagesc(Bx);
subplot(2,3,2+1*3); imagesc(abs(Bk));
subplot(2,3,3+1*3); imagesc(angle(Bk));
end;%if (disp_flag);



