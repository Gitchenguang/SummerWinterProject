function [Xpset,Ypset]=pointsset(step, Axy,Bxy,Cxy,Dxy,Exy,Fxy  )
% [Xpset,Ypset]=pointsset( Axy,Bxy,Cxy,Dxy,Exy,Fxy  )

AB=Bxy-Axy;BC=Cxy-Bxy;CD=Dxy-Cxy;AF=Fxy-Axy;FE=Exy-Fxy;ED=Dxy-Exy;
kAB=AB(1,2)/AB(1,1);kAF=AF(1,2)/AF(1,1);kED=ED(1,2)/ED(1,1);kCD=CD(1,2)/CD(1,1);

length=abs(FE(1,1));

index=1;

for i=step/2:step:2*length-step/2
    for j=step/2:step:2*length-step/2
        % 对i分区间讨论[0,25]
        if i<=25
            if (j>( kAB*(i-Axy(1,1))+Axy(1,2)))&&(j<( kAF*(i-Axy(1,1))+Axy(1,2)) )
                Xpset(index)=i;
                Ypset(index)=j;
                index=index+1;
            end
        elseif i<=75
            if j<=Fxy(1,2)
                Xpset(index)=i;
                Ypset(index)=j;
                index=index+1;
            end
        elseif i<=100
            if (j>( kCD*(i-Cxy(1,1))+Cxy(1,2) ))&&(j<( kED*(i-Exy(1,1))+Exy(1,2)) )
                Xpset(index)=i;
                Ypset(index)=j;
                index=index+1;
            end
        end
    end
end
