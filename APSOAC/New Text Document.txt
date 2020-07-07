clc;clear all;

f2=ones(18,9).*inf;
f2min=ones(18,1).*inf;
P1min=zeros(18,1);
P2min=zeros(18,1);
D12=zeros(18,1);
Ans12=zeros(18,1);

f3=ones(27,9).*inf;
f3min=ones(27,1).*inf;
P3min=zeros(27,1);
D123=zeros(27,1);
Ans123=zeros(27,1);

Powers=[0;50;75;100;125;150;175;200;225];
Costs=[ inf    inf  inf;
        810    750  806;
        1355   1155 1108.5;
        1460   1360 1411;
        1772.5 1655 1704.5;
        2085   1950 1998;
        2427.5 inf  2358;
        2760   inf  inf;
        inf    inf  inf];

i1min=2;i2min=2;i3min=2;
i1max=8;i2max=6;i3max=7;
D12min=100;D12max=350;
D123min=300;D123max=325;

for c=i2min:i2max
    for r=i1min:i1max
        sumP=Powers(r)+Powers(c);
        if (sumP>=D12min)&&(sumP<=D12max)
            f2(r+c,c)=Costs(r,1)+Costs(c,2);
            if (f2(r+c,c)==min(f2(r+c,:)))&&(f2(r+c,c)~=inf)
                f2min(r+c,1)=min(f2(r+c,:));
                P2min(r+c,1)=Powers(c,1);
                P1min(r+c,1)=Powers(r,1);
                D12(r+c,1)=Powers(r)+Powers(c);
                Ans12(r+c,1)=D12(r+c,1);
                Ans12(r+c,2)=f2min(r+c,1);
                Ans12(r+c,3)=P2min(r+c,1);
                Ans12(r+c,4)=P1min(r+c,1);
            end
        end
    end
end

for c=i3min:i3max
    for r=1:length(f2min)
        sumP=D12(r)+Powers(c);
        if (sumP>=D123min)&&(sumP<=D123max)
            f3(r+c,c)=f2min(r,1)+Costs(c,3);
            if (f3(r+c,c)==min(f3(r+c,:)))&&(f3(r+c,c)~=inf)
                f3min(r+c,1)=min(f3(r+c,:));
                P3min(r+c,1)=Powers(c,1);
                D123(r+c,1)=D12(r)+Powers(c,1);
                Ans123(r+c,1)=D123(r+c,1);
                Ans123(r+c,2)=f3min(r+c,1);
                Ans123(r+c,3)=P3min(r+c,1);
                Ans123(r+c,4)=P2min(r,1);
                Ans123(r+c,5)=P1min(r,1);
            end
        end
    end
end

Ans12
Ans123