function [ amp ] = mag3( a,b,c )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[m,n]=size(a);
for ro=1:m
    for co=1:n
        amp(ro,co)=sqrt(a(ro,co)*a(ro,co)+b(ro,co)*b(ro,co)+c(ro,co)*c(ro,co));
    end
end

