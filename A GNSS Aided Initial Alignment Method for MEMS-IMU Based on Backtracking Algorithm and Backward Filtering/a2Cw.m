function [ Cw ] = a2Cw( a )
%   a2Cw ŷ������תCw
    Cw = [cos(a(2)) 0 -sin(a(2))*cos(a(1))
              0     1    sin(a(1))        
          sin(a(2)) 0  cos(a(2))*cos(a(1))  ];
end

