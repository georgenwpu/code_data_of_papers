function [ Cn_n ] = a2Cn_n( a )
%a2Cn_n Å·À­Îó²î½Ç×ª×ËÌ¬Îó²îÕó
    ax = a(1);  sax = sin(ax);  cax = cos(ax);
    ay = a(2);  say = sin(ay);  cay = cos(ay);
    az = a(3);  saz = sin(az);  caz = cos(az);
    Cn_n = [cay*caz-say*sax*saz  cay*saz+say*sax*caz -say*cax
               -cax*saz               cax*caz         sax
           say*caz+cay*sax*saz  say*saz-cay*sax*caz cay*cax];
end

