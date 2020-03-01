clc;clear all;

%%
Ngen=3;
Pmax=[600 400 200];
Pmin=[150 100 050];
Pload=850;
divisions=50;

%%
for i=1:Ngen
    range(i)=Pmax(i)-Pmin(i);
    dP(i)=range(i)/divisions;
end

%%
for i=1:Ngen
    for k=1:divisions
        sPmin(i,k)=Pmin(i)+(k-1)*dP(i);
        sPmax(i,k)=Pmin(i)+k    *dP(i);
        sFmin(i,k)=F(i,sPmin(i,k));
        sFmax(i,k)=F(i,sPmax(i,k));
        s    (i,k)=(sFmax(i,k)-sFmin(i,k))/dP(i);
    end
end
ordered_s=sort(transpose(reshape(s,[],1)));

%%
Pgen=[Pmin(1) Pmin(2) Pmin(3)];
eps=abs(sum(Pgen)-Pload);
threshold=1;
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
    plot([50,600],[lamda,lamda])
    %plot(Pgen,lamda,'o')
    
    for rep=1:2
        for i=1:Ngen
            for k=1:divisions
                
                if (s(i,k)<lamda)&&(max(s(i,1:k))==s(i,k))
                    Pgen(i)=sPmax(i,k);
                end
                if (s(i,k)>lamda)&&(min(s(i,k:divisions))==s(i,k))
                    Pgen(i)=sPmin(i,k);
                end
                
                for b=1:Ngen
                    if (Pgen(b)<Pmin(b))
                        Pgen(b)=Pmin(b);
                    end
                    if (Pgen(b)>Pmax(b))
                        Pgen(b)=Pmax(b);
                    end
                    
                    for p=1:divisions
                        if (s(b,p)==lamda)
                            Pgen;
                            t=Pload-sum(Pgen)+Pgen(b);
                            if ((t<=sPmax(b,p))&&(t>=sPmin(b,p)))
                                Pgen(b)=t;
                            end
                        end
                    end
                end
                Pgen;
                eps=abs(sum(Pgen)-Pload);
                eps_his(n*i*k)=eps;
                TotalCost=0;
                for m=1:Ngen
                    TotalCost=TotalCost+F(m,Pgen(m));
                end
                
                if(eps<threshold)&&(TotalCost<OldCostMin)
                    Ans=Pgen;
                    TotalCost;
                    Oldeps=eps;
                    OldCostMin=TotalCost;
                end
            end
        end
    end
    
    plot(Pgen,lamda,'o')
    Pgen=[0 0 0];
    %     Pgen=[Pmin(1) Pmin(2) Pmin(3)];
end

Pgen=Ans
TotalCost=0;
for i=1:Ngen
    TotalCost=TotalCost+F(i,Pgen(i));
end
TotalCost
