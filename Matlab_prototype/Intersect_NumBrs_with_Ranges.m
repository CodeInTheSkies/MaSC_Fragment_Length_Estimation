% This function accepts a set of sorted numbers and a set of ranges (or intervals)
% and outputs only the numbers that fall within the ranges (intervals)

function OtptNos = Intersect_NumBrs_with_Ranges(InptNos, RngStrts, RngEnds)

Inptcnt = 1;  Rngcnt = 1;  % Assigning separate variables for counting

% Preallocating the output variable so as to make the execution faster
OtptNos = zeros(1,length(InptNos));  Otptcnt = 0;  

while ~(Inptcnt > length(InptNos)) && ~(Rngcnt > length(RngStrts))
%     
if InptNos(Inptcnt) < RngStrts(Rngcnt)
    Inptcnt = Inptcnt + 1;
    continue;
end
% 
if (InptNos(Inptcnt) >= RngStrts(Rngcnt)) && (InptNos(Inptcnt) <= RngEnds(Rngcnt))
    Otptcnt = Otptcnt + 1;
    OtptNos(Otptcnt) = InptNos(Inptcnt);
    Inptcnt = Inptcnt + 1;
    continue;
end
% 
if InptNos(Inptcnt) > RngEnds(Rngcnt)
    Rngcnt = Rngcnt + 1;
    continue;
end
% 
end
% Reassigning the output variable to trim any training zeros from preallocation
OtptNos = OtptNos(1:Otptcnt);
