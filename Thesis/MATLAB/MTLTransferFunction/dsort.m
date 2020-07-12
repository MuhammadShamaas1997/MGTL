function [ c,d ] = dsort( a,b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for in=1:length(a)
    for in1=1:length(a)-1
        if(a(in1)>a(in1+1))
            t=a(in1);
            c(in1)=a(in1+1);
            a(in1)=a(in1+1);
            c(in1+1)=t;
            a(in1+1)=t;
            
            t=b(in1);
            d(in1)=b(in1+1);
            b(in1)=b(in1+1);
            d(in1+1)=t;
            b(in1+1)=t;
        end
    end
end


end

