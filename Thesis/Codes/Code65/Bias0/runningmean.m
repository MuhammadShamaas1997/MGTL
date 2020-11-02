function [ out ] = runningmean( inm,num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
extra=((num-1)/2);
out=zeros(1,length(inm));
inm=[zeros(1,extra) inm zeros(1,extra)];
for pl=(1+extra):(length(inm)-extra)
sumn=mean(inm((pl-extra):(pl+extra)));
    out(pl-extra)=sumn/num;
end

end

