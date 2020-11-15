function h = activation_simple(type, val, constraint, ths)

%Majority constraint
if (type == '>')
    if (val >= constraint+ths)
        h = 0;
    elseif (val < constraint+ths && val > constraint)
        a = (val-constraint)/ths;
        h = .5*(cos(a*pi)+1);%Se usiamo questa: 6*a^5-15*a^4+10*a^3, attivazione: a = (constraint+ths-val)/ths;
    else
        h = 1;
    end
end

if (type == '<')
    if (val < constraint-ths)
        h = 0;
    else
        if (val <= constraint && val >= constraint-ths)
            a = (val-constraint)/ths;
            h = .5*(cos(a*pi)+1);%Se usiamo questa: 6*a^5-15*a^4+10*a^3, attivazione: a = (constraint+ths-val)/ths;
        else
            h = 1;
        end
    end
end

