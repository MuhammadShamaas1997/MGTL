function [ c ] = rs( i,p )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (i==1)
    if (p<=600)&&(p>=150)
        c=7.92+2*0.001562*p;
    else
        c=0;
    end
elseif (i==2)
    if (p<=400)&&(p>=100)
    c=7.85+2*0.00194*p;
    else
        c=0;    
    end
elseif (i==3)
    if (p<=200)&&(p>=50)
    c=7.97+2*0.00482*p;
    else
        c=0;
    end
end

end

