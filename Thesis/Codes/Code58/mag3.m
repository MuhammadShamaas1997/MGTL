function [ amp ] = mag3( a,b,c )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(a)
    amp(i)=sqrt(a(i)*a(i)+b(i)*b(i)+c(i)*c(i));
end
end

