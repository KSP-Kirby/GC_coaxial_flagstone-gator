%load('flow.mat')
front = 0;
if(front == 1 )
    constrast1 = 1;
    offset1 = 0;
    constrast2 = 1.5;
    offset2 = 0;
    constrast3 = 1.5;
    offset3 = 0;
    imgs = cast(uvi_f{1},'double')/256;
    imgs1 = imgs(:,:,1);
    imgs2 = imgs1;
    imgs2(:,:,2) = (imgs(:,:,2)-offset2)*constrast2;
    imgs2(:,:,3) = (imgs(:,:,3)-offset3)*constrast3;
else
    constrast1 = 1;
    offset1 = 0;
    constrast2 = 1.5;
    offset2 = 0;
    constrast3 = 1.5;
    offset3 = 0;
    imgs = cast(uvi_b{1},'double')/256;
    imgs1 = imgs(:,:,1);
    imgs2 = imgs1;
    imgs2(:,:,2) = (imgs(:,:,2)-offset2)*constrast2;
    imgs2(:,:,3) = (imgs(:,:,3)-offset3)*constrast3;    
end

imtool(mirrorHorz(imgs2))