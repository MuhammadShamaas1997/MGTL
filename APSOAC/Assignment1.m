clc;clear all;
Ngen=3;
Pmax=[600 400 200];
Pmin=[150 100 050];
Pload=850;
divisions=50;

for i=1:Ngen
    range(i)=Pmax(i)-Pmin(i);
    dP(i)=range(i)/divisions;
end

for i=1:Ngen
    for k=1:divisions
        sPmin(i,k)=Pmin(i)+(k-1)*dP(i);
        sPmax(i,k)=Pmin(i)+k    *dP(i);
        sFmin(i,k)=F(i,sPmin(i,k));
        sFmax(i,k)=F(i,sPmax(i,k));
        s    (i,k)=(sFmax(i,k)-sFmin(i,k))/dP(i);
        s  (i,k+1)=(sFmax(i,k)-sFmin(i,k))/dP(i);
        
    end
end

ordered_s=sort(transpose(reshape(s,[],1)));


Pgen=[Pmin(1) Pmin(2) Pmin(3)];
eps=abs(sum(Pgen)-Pload);
threshold=2;
iter=1;
maxiter=Ngen*divisions;
TotalCost=0;
found=0;
Ans=[];
OldCostMin=1e5;
Oldeps=Pload;

for n=1:maxiter
    lamda=ordered_s(n);
    hold on
%     stairs(sPmin(1,:),s(1,:),'r')
%     stairs(sPmin(2,:),s(2,:),'b')
%     stairs(sPmin(3,:),s(3,:),'g')
    plot([50,600],[lamda,lamda],'*')
    plot(Pgen,lamda,'o')
    
             
    for i=1:Ngen
        for k=1:divisions
            if (min(s(i,:))>lamda)
                Pgen(i)=Pmin(i);
            end
            if (max(s(i,:))<lamda)
                Pgen(i)=Pmax(i);
            end
            if (s(i,k)<lamda)&&(s(i,k+1)>lamda)
                Pgen(i)=sPmax(i,k);
            end
            if (s(i,k)==lamda)
                t=Pload-sum(Pgen)+Pgen(i);
                Pgen(i)=sPmin(i,k);
                if ((t<=sPmax(i,k))&&(t>=sPmin(i,k)))
                Pgen(i)=t;
                end
            end
            %Pgen
            eps=abs(sum(Pgen)-Pload);
            eps_his(iter)=eps;
            TotalCost=0;
            for m=1:Ngen
                TotalCost=TotalCost+F(m,Pgen(m));
            end
            
             if(eps<Oldeps)
                 Ans=Pgen;
                 Oldeps=eps;
                 OldCostMin=TotalCost;
             end
            iter=iter+1;
        end
    end
end

Pgen=Ans
TotalCost=0;
for i=1:Ngen
    TotalCost=TotalCost+F(i,Pgen(i));
end
TotalCost
