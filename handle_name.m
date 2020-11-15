function RTPflag = handle_name(name)
%HANDLE_NAME Summary of this function goes here
%   Detailed explanation goes here

RTPflag = 0;
if (strcmp(name,"iCAT-SR"))
    RTPflag = 5;
else
    if (strcmp(name,"RP"))
        RTPflag = 1;
    else
        if (strcmp(name,"iCAT-TPC"))
            RTPflag = 2;
        else
            if (strcmp(name,"iCAT-RP"))
                RTPflag = 3;
            else
                if (strcmp(name,"ST"))
                    RTPflag = 4;
                else
                    if (strcmp(name,"iCAT-conservative-RP"))
                        RTPflag = 6;
                    else
                        disp("Default SR");
                    end
                end
            end
        end
    end
end

