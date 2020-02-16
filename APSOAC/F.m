function [ c ] = F( i,p )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (i==1)
    if (p<=600)&&(p>=150)
        c=561+7.92*p+0.001562*p*p;
    else
        c=0;
    end
elseif (i==2)
    if (p<=400)&&(p>=100)
    c=310+7.85*p+0.00194*p*p;
    else
        c=0;    
    end
elseif (i==3)
    if (p<=200)&&(p>=50)
    c=78+7.97*p+0.00482*p*p;
    else
        c=0;
    end
end

end

