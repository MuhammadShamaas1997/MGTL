function [ c ] = ph( a,b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(a)
    c(i)=atan(b(i)/a(i));
end

end

