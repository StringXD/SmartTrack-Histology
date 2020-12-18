% implanted tracks rotation video
fh = open('demo2.fig');
WriterObj = VideoWriter('allenCCFsliceMovie');
WriterObj.FrameRate=30;
open(WriterObj);
nFr = 720;
% design a view trajectory
viewPoint = zeros(nFr,2);
viewPoint(1:360,1) = (1:360)';
viewPoint(361:540,1) = (1:180)';
viewPoint(541:720,1) = (179:-1:0)';
viewPoint(361:450,2) = (1:90)';
viewPoint(451:630,2) = (89:-1:-90)';
viewPoint(631:720,2) = (-89:0)';
for i = 1:nFr
    view(viewPoint(i,1),viewPoint(i,2));
    frame = getframe(fh);
    writeVideo(WriterObj,frame);
end
close(WriterObj);