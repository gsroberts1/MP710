fullfig; imshow(stiff,[]);
ROI1 = drawrectangle('Color','g','Label','ROI1');
title('True Stiffness Map');
i1 = floor(ROI1.Position);          % Position of the ROI [xmin, ymin, width, height]
                                    % xmin and ymin = upper left corner of ROI
rows = i1(2):(i1(2)+i1(4));
cols = i1(1):(i1(1)+i1(3));
save rows.mat rows
save cols.mat cols