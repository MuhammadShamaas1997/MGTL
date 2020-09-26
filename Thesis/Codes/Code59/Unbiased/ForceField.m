clc;clear all;
f=fopen('bx-000050.25.h5');

ex=h5read('ex-000050.25.h5','/ex.r');
ey=h5read('ey-000050.25.h5','/ey.r');
ez=h5read('ez-000050.25.h5','/ez.r');
hold on;
for i=1:length(ex)
    for j=1:33
              ex(i,j)=ex(i,j)/1e5;
              ez(i,j)=ez(i,j)/1e5;
%             e(i,j)=ex(i,j)*ex(i,j)+ey(i,j)*ey(i,j)+ez(i,j)*ez(i,j);
            plot([i-ex(i,j),i+ex(i,j)],[j-ez(i,j),j+ez(i,j)]);
            plot([i+ex(i,j),i+ex(i,j)],[j+.5*ez(i,j),j+ez(i,j)]);
            plot([i+.5*ex(i,j),i+ex(i,j)],[j+ez(i,j),j+ez(i,j)]);
    end
end

   
dx=h5read('dx-000050.25.h5','/dx.r');
dy=h5read('dy-000050.25.h5','/dy.r');
dz=h5read('dz-000050.25.h5','/dz.r');
hold on;
for i=1:length(dx)
    for j=1:length(dy)
%             d(i,j)=dx(i,j)*dx(i,j)+dy(i,j)*dy(i,j)+dz(i,j)*dz(i,j);
%             plot([i-.50*dx(i,j),i+.50*dx(i,j)],[j-.50*dy(i,j),j+.50*dy(i,j)]);
%             plot([i+.50*dx(i,j),i+.50*dx(i,j)],[j+.25*dy(i,j),j+.50*dy(i,j)]);
%             plot([i+.25*dx(i,j),i+.50*dx(i,j)],[j+.50*dy(i,j),j+.50*dy(i,j)]);
    end
end


hx=h5read('hx-000050.25.h5','/hx.r');
hy=h5read('hy-000050.25.h5','/hy.r');
hz=h5read('hz-000050.25.h5','/hz.r');
hold on;
for i=1:length(hx)
    for j=1:length(hy)
%             h(i,j)=hx(i,j)*hx(i,j)+hy(i,j)*hy(i,j)+hz(i,j)*hz(i,j);
%             plot([i-.50*hx(i,j),i+.50*hx(i,j)],[j-.50*hy(i,j),j+.50*hy(i,j)]);
%             plot([i+.50*hx(i,j),i+.50*hx(i,j)],[j+.25*hy(i,j),j+.50*hy(i,j)]);
%             plot([i+.25*hx(i,j),i+.50*hx(i,j)],[j+.50*hy(i,j),j+.50*hy(i,j)]);
    end
end


bx=h5read('bx-000050.25.h5','/bx.r');
by=h5read('by-000050.25.h5','/by.r');
bz=h5read('bz-000050.25.h5','/bz.r');
hold on;
for i=1:length(by)
    for j=1:length(bz)
%             b(i,j)=bx(i,j)*bx(i,j)+by(i,j)*by(i,j)+bz(i,j)*bz(i,j);
            %bx(i,j)=bx(i,j)/max(bx);
%             by(i,j)=by(i,j)/1e5;
%             bz(i,j)=bz(i,j)/1e5;
%             plot([i-by(i,j),i+by(i,j)],[j-bz(i,j),j+bz(i,j)]);
%             plot([i+by(i,j),i+by(i,j)],[j+0.5*bz(i,j),j+bz(i,j)]);
%             plot([i+0.5*by(i,j),i+by(i,j)],[j+bz(i,j),j+bz(i,j)]);
    end
end


sx=h5read('sx-000050.25.h5','/sx');
sy=h5read('sy-000050.25.h5','/sy');
sz=h5read('sz-000050.25.h5','/sz');

hold on;
for i=1:length(sy)
    for j=1:length(sz)
%             s(i,j)=sx(i,j)*sx(i,j)+sy(i,j)*sy(i,j)+sz(i,j)*sz(i,j);
%             plot([i-.0000050*sy(i,j),i+.0000050*sy(i,j)],[j-.0000050*sz(i,j),j+.0000050*sz(i,j)]);
%             plot([i+.0000050*sy(i,j),i+.0000050*sy(i,j)],[j+.0000025*sz(i,j),j+.0000050*sz(i,j)]);
%             plot([i+.0000025*sy(i,j),i+.0000050*sy(i,j)],[j+.0000050*sz(i,j),j+.0000050*sz(i,j)]);
    end
end