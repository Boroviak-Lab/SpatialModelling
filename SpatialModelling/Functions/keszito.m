function [ x_ad,y_ad,z_ad, kpi ] = keszito( f )
%KESZITO Summary of this function goes here
szam=0;
x_ad=zeros(3000,1);
y_ad=zeros(3000,1);
z_ad=zeros(3000,1);
kpi=zeros(f.pieces,1);
for i=1:f.pieces %Loop over pieces
    kpi(i)=szam+1; %Counter
    for t=0:0.1:abs(f.breaks(1,i)-f.breaks(1,i+1)) %Break current break into increments of .1?
        szam=szam+1;
        
        %Explicit order 4?
        x_ad(szam,:)=f.coefs(3*i-2,1)*t^3 + f.coefs(3*i-2,2)*t^2 + f.coefs(3*i-2,3)*t^1 + f.coefs(3*i-2,4)*t^0;
        
        y_ad(szam,:)=f.coefs(3*i-1,1)*t^3+f.coefs(3*i-1,2)*t^2+f.coefs(3*i-1,3)*t^1+f.coefs(3*i-1,4)*t^0;

        z_ad(szam,:)=f.coefs(3*i,1)*t^3+f.coefs(3*i,2)*t^2+f.coefs(3*i,3)*t^1+f.coefs(3*i,4)*t^0;
        
    end
    if i==f.pieces && t~=abs(f.breaks(1,i)-f.breaks(1,i+1)) % this part is for the last point (curve length / t is not an integer usually)
        t=abs(f.breaks(1,i)-f.breaks(1,i+1));
        szam=szam+1;
        x_ad(szam,:)=f.coefs(3*i-2,1)*t^3+f.coefs(3*i-2,2)*t^2+f.coefs(3*i-2,3)*t^1+f.coefs(3*i-2,4)*t^0;
        y_ad(szam,:)=f.coefs(3*i-1,1)*t^3+f.coefs(3*i-1,2)*t^2+f.coefs(3*i-1,3)*t^1+f.coefs(3*i-1,4)*t^0;
        z_ad(szam,:)=f.coefs(3*i,1)*t^3+f.coefs(3*i,2)*t^2+f.coefs(3*i,3)*t^1+f.coefs(3*i,4)*t^0;

    end
end
kpi(i+1)=szam;
x_ad=x_ad(1:szam,1);
y_ad=y_ad(1:szam,1); 
z_ad=z_ad(1:szam,1); 

end