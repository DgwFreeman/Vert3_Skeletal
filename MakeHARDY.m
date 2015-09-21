function HARDY = MakeHARDY()

%% HARDY will be binding sites or actin node numbers:
%HARDY=zeros(ThinRow_Index, floor(NumSites/3) );
HARDY=[
% For actin node a)
2:3:91;
3:3:91;
4:3:91;
% For actin node b)
93:3:182;
94:3:182;
95:3:182;
% For actin node c)
184:3:273;
185:3:273;
186:3:273;
% For actin node d)
275:3:364;
276:3:364;
 277:3:364;
% For actin node e)
366:3:455;
367:3:455;
368:3:455;
% For actin node f)
457:3:546;
458:3:546;
459:3:546;
% For actin node g)
548:3:637;
549:3:637;
550:3:637;
% For actin node h)
639:3:728;
640:3:728;
641:3:728];