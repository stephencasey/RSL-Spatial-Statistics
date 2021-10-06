function a1String = idx2A1(idx)

alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';

if idx < 27
    a1String = alphabet(idx);
else
    idx2=rem(idx,26);
    if idx2 == 0
        a1String=[alphabet(floor(idx/26)-1),'Z'];
    else
        a1String=[alphabet(floor(idx/26)),alphabet(idx2)];
    end
end