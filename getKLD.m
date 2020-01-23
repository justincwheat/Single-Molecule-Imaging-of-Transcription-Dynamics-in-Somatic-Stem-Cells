function y = getKLD(Obsdist, Simdist)

    
    y= sum(nansum(Obsdist.*log2(Simdist./Obsdist)));
    
   