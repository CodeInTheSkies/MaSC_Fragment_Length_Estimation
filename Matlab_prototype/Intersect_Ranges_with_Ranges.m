% This function accepts two sets of sorted ranges (or intervals)
% and outputs total number of overlaps between them (intersections)

function OvrLaps = Intersect_Ranges_with_Ranges(Rng1Strts, Rng1Ends, Rng2Strts, Rng2Ends)

Rng1cnt = 1;  Rng2cnt = 1;  % Assigning separate variables for counting

OvrLaps = 0; 

while ~(Rng1cnt > length(Rng1Strts)) && ~(Rng2cnt > length(Rng2Strts))
%     
if Rng1Ends(Rng1cnt) < Rng2Strts(Rng2cnt)
    Rng1cnt = Rng1cnt + 1;
elseif Rng1Strts(Rng1cnt) > Rng2Ends(Rng2cnt)
    Rng2cnt = Rng2cnt + 1;
else
    OvrLaps = OvrLaps + min(Rng1Ends(Rng1cnt), Rng2Ends(Rng2cnt)) - max(Rng1Strts(Rng1cnt), Rng2Strts(Rng2cnt)) + 1;
    if Rng1Ends(Rng1cnt) < Rng2Ends(Rng2cnt)
        Rng1cnt = Rng1cnt + 1;
    elseif Rng1Ends(Rng1cnt) > Rng2Ends(Rng2cnt)
        Rng2cnt = Rng2cnt + 1;
    else
        Rng1cnt = Rng1cnt + 1;
        Rng2cnt = Rng2cnt + 1;
    end        
end
%
end