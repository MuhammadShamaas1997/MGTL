clc;clear all;
f=fopen('bx-001000.00.h5');

ex=h5read('ex-001000.00.h5','/ex.r');
ey=h5read('ey-001000.00.h5','/ey.r');
ez=h5read('ez-001000.00.h5','/ez.r');
hold on;
for i=1:length(ex)
    for j=1:length(ey)
            e(i,j)=ex(i,j)*ex(i,j)+ey(i,j)*ey(i,j)+ez(i,j)*ez(i,j);
%             plot([i-50*ex(i,j),i+50*ex(i,j)],[j-50*ey(i,j),j+50*ey(i,j)]);
%             plot([i+50*ex(i,j),i+50*ex(i,j)],[j+25*ey(i,j),j+50*ey(i,j)]);
%             plot([i+25*ex(i,j),i+50*ex(i,j)],[j+50*ey(i,j),j+50*ey(i,j)]);
    end
end

   
dx=h5read('dx-001000.00.h5','/dx.r');
dy=h5read('dy-001000.00.h5','/dy.r');
dz=h5read('dz-001000.00.h5','/dz.r');
hold on;
for i=1:length(dx)
    for j=1:length(dy)
            d(i,j)=dx(i,j)*dx(i,j)+dy(i,j)*dy(i,j)+dz(i,j)*dz(i,j);
%             plot([i-.50*dx(i,j),i+.50*dx(i,j)],[j-.50*dy(i,j),j+.50*dy(i,j)]);
%             plot([i+.50*dx(i,j),i+.50*dx(i,j)],[j+.25*dy(i,j),j+.50*dy(i,j)]);
%             plot([i+.25*dx(i,j),i+.50*dx(i,j)],[j+.50*dy(i,j),j+.50*dy(i,j)]);
    end
end


hx=h5read('hx-001000.00.h5','/hx.r');
hy=h5read('hy-001000.00.h5','/hy.r');
hz=h5read('hz-001000.00.h5','/hz.r');
hold on;
for i=1:length(hx)
    for j=1:length(hy)
            h(i,j)=hx(i,j)*hx(i,j)+hy(i,j)*hy(i,j)+hz(i,j)*hz(i,j);
%             plot([i-.50*hx(i,j),i+.50*hx(i,j)],[j-.50*hy(i,j),j+.50*hy(i,j)]);
%             plot([i+.50*hx(i,j),i+.50*hx(i,j)],[j+.25*hy(i,j),j+.50*hy(i,j)]);
%             plot([i+.25*hx(i,j),i+.50*hx(i,j)],[j+.50*hy(i,j),j+.50*hy(i,j)]);
    end
end


bx=h5read('bx-001000.00.h5','/bx.r');
by=h5read('by-001000.00.h5','/by.r');
bz=h5read('bz-001000.00.h5','/bz.r');
hold on;
for i=1:length(bx)
    for j=1:length(by)
            b(i,j)=bx(i,j)*bx(i,j)+by(i,j)*by(i,j)+bz(i,j)*bz(i,j);
%             plot([i-.50*bx(i,j),i+.50*bx(i,j)],[j-.50*by(i,j),j+.50*by(i,j)]);
%             plot([i+.50*bx(i,j),i+.50*bx(i,j)],[j+.25*by(i,j),j+.50*by(i,j)]);
%             plot([i+.25*bx(i,j),i+.50*bx(i,j)],[j+.50*by(i,j),j+.50*by(i,j)]);
    end
end


sx=h5read('sx-001000.00.h5','/sx');
sy=h5read('sy-001000.00.h5','/sy');
sz=h5read('sz-001000.00.h5','/sz');
hold on;
for i=1:length(sx)
    for j=1:length(sy)
            s(i,j)=sx(i,j)*sx(i,j)+sy(i,j)*sy(i,j)+sz(i,j)*sz(i,j);
            plot([i-.050*sx(i,j),i+.050*sx(i,j)],[j-.050*sy(i,j),j+.050*sy(i,j)]);
            plot([i+.050*sx(i,j),i+.050*sx(i,j)],[j+.025*sy(i,j),j+.050*sy(i,j)]);
            plot([i+.025*sx(i,j),i+.050*sx(i,j)],[j+.050*sy(i,j),j+.050*sy(i,j)]);
    end
end