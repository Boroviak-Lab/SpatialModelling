function [Pr] = Quaternion3(theta,u,x);

% Quaternion rotation.
%    Copyright (C) 2010  Christopher Penfold

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

si      = length(x(:,1));
vv      = u*sin(theta/2);
v       = vv(ones(si,1),:);
s       = cos(theta/2) + zeros(si,1);
qP      = [ s(:,ones(1,3)).*x +  ...
             [v(:,2).*x(:,3) - v(:,3).*x(:,2),v(:,3).*x(:,1) - v(:,1).*x(:,3),v(:,1).*x(:,2) - v(:,2).*x(:,1)]];
       In1 = - sum(conj(v).*x,2);
crossing2 = [qP(:,2).*(-v(:,3)) - qP(:,3).*(-v(:,2)),qP(:,3).*(-v(:,1)) - qP(:,1).*(-v(:,3)),qP(:,1).*(-v(:,2)) - qP(:,2).*(-v(:,1))];
        Pr      = [In1(:,ones(1,3)).*(-v)+ s(:,ones(1,3)).*qP + crossing2];
