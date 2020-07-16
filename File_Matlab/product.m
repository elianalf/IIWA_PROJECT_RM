function t=product(a,b)
x=a(2)*b(3)-a(3)*b(2);
y=a(3)*b(1)-a(1)*b(3);
z=a(1)*b(2)-a(2)*b(1);
t=[x,y,z]';
end