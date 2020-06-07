clc;clear all;
f=fopen('bx-000200.50.h5');

ex=h5read('ex-000200.50.h5','/ex.r');
ey=h5read('ey-000200.50.h5','/ey.r');
ez=h5read('ez-000200.50.h5','/ez.r');
hold on;
for i=1:length(ey)
    for j=1:length(ez)
            ex(i,j)=1000*ex(i,j);
            ey(i,j)=1000*ey(i,j);
            ez(i,j)=1000*ez(i,j);            
    end
end

for i=20:2:length(ey)-20
    for j=20:2:length(ez)-20
%             plot([i-ey(i,j),i+ey(i,j)],[j-ez(i,j),j+ez(i,j)]);
%             plot([i+ey(i,j),i+ey(i,j)],[j+.5*ez(i,j),j+ez(i,j)]);
%             plot([i+.5*ey(i,j),i+ey(i,j)],[j+ez(i,j),j+ez(i,j)]);
    end
end

   
dx=h5read('dx-000200.50.h5','/dx.r');
dy=h5read('dy-000200.50.h5','/dy.r');
dz=h5read('dz-000200.50.h5','/dz.r');
hold on;
for i=1:length(dy)
    for j=1:length(dz)
%             d(i,j)=dx(i,j)*dx(i,j)+dy(i,j)*dy(i,j)+dz(i,j)*dz(i,j);
            dx(i,j)=1000*dx(i,j);
            dy(i,j)=1000*dy(i,j);
            dz(i,j)=1000*dz(i,j);            

    end
end

for i=20:5:length(dy)-20
    for j=20:5:length(dz)-20
            plot([i-dy(i,j),i+dy(i,j)],[j-dz(i,j),j+dz(i,j)]);
            plot([i+dy(i,j),i+dy(i,j)],[j+.5*dz(i,j),j+dz(i,j)]);
            plot([i+.5*dy(i,j),i+dy(i,j)],[j+dz(i,j),j+dz(i,j)]);
    end
end

hx=h5read('hx-000200.50.h5','/hx.r');
hy=h5read('hy-000200.50.h5','/hy.r');
hz=h5read('hz-000200.50.h5','/hz.r');
hold on;
for i=1:length(hx)
    for j=1:length(hy)
            h(i,j)=hx(i,j)*hx(i,j)+hy(i,j)*hy(i,j)+hz(i,j)*hz(i,j);
            hx(i,j)=4000*hx(i,j);
            hy(i,j)=4000*hy(i,j);
            hz(i,j)=4000*hz(i,j);            
    end
end

for i=20:5:length(hx)-20
    for j=20:5:length(hy)-20            
%             plot([i-hy(i,j),i+hy(i,j)],[j-hz(i,j),j+hz(i,j)]);
%             plot([i+hy(i,j),i+hy(i,j)],[j+.5*hz(i,j),j+hz(i,j)]);
%             plot([i+.5*hy(i,j),i+hy(i,j)],[j+hz(i,j),j+hz(i,j)]);
    end
end


bx=h5read('bx-000200.50.h5','/bx.r');
by=h5read('by-000200.50.h5','/by.r');
bz=h5read('bz-000200.50.h5','/bz.r');
hold on;
for i=1:length(by)
    for j=1:length(bz)
            b(i,j)=bx(i,j)*bx(i,j)+by(i,j)*by(i,j)+bz(i,j)*bz(i,j);
%             bx(i,j)=100*bx(i,j)/max(max(bx));
%             by(i,j)=100*by(i,j)/max(max(by));
%             bz(i,j)=100*bz(i,j)/max(max(bz));            
            bx(i,j)=100*bx(i,j);
            by(i,j)=100*by(i,j);
            bz(i,j)=100*bz(i,j);            
    end
end

for i=1:5:length(by)
    for j=1:5:length(bz)
%             plot([i-by(i,j),i+by(i,j)],[j-bz(i,j),j+bz(i,j)]);
%             plot([i+by(i,j),i+by(i,j)],[j+0.5*bz(i,j),j+bz(i,j)]);
%             plot([i+0.5*by(i,j),i+by(i,j)],[j+bz(i,j),j+bz(i,j)]);
    end
end

sx=h5read('sx-000200.50.h5','/sx');
sy=h5read('sy-000200.50.h5','/sy');
sz=h5read('sz-000200.50.h5','/sz');

hold on;
for i=1:length(sy)
    for j=1:length(sz)
%             s(i,j)=sx(i,j)*sx(i,j)+sy(i,j)*sy(i,j)+sz(i,j)*sz(i,j);
%             plot([i-.0000050*sy(i,j),i+.0000050*sy(i,j)],[j-.0000050*sz(i,j),j+.0000050*sz(i,j)]);
%             plot([i+.0000050*sy(i,j),i+.0000050*sy(i,j)],[j+.0000025*sz(i,j),j+.0000050*sz(i,j)]);
%             plot([i+.0000025*sy(i,j),i+.0000050*sy(i,j)],[j+.0000050*sz(i,j),j+.0000050*sz(i,j)]);
    end
end