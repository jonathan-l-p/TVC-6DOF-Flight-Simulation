clc;

% Plot rotation

plot = true;

% compile pry and thrust into one vector

if ((length(out.p.Time) == length(out.r.Time)) && (length(out.r.Time) == length(out.y.Time)))
    r_time = out.r.Time;
    thrust_time = out.T.Time;

    match = []; % vector of the coordinates of thrust time that match r time
    
    for i = 1:length(r_time)
        for j = 1:length(thrust_time)
            
            if thrust_time(j) == r_time(i)
                match(i) = j;
            end
            
        end
    end
    
    thrust2_time = out.T.Time(match,:);
    thrust2 = out.T.Data(match,:);
    
    
    % plot
    if sum(thrust2_time == r_time) ~= length(r_time)
        disp('Error: vector size')
    elseif plot == true
        
        figure(1)
        
        for k = 1:70 % k = 1:length(r_time)
            p0 = [0,0,0];
            p1 = [out.p.Data(k,1), out.p.Data(k,3), out.p.Data(k,2)];
            vectarrow(p0,p1,'k')
            hold on
            r0 = [0,0,0];
            r1 = 2*[out.r.Data(k,1), out.r.Data(k,3), out.r.Data(k,2)];
            vectarrow(r0,r1,'k')
            hold on
            y0 = [0,0,0];
            y1 = [out.y.Data(k,1), out.y.Data(k,3), out.y.Data(k,2)];
            vectarrow(y0,y1,'k')
            hold on
            thrust0 = [0,0,0];
            thrust1 = -[thrust2(k,1), thrust2(k,3), thrust2(k,2)]/30;
            vectarrow(thrust0,thrust1,'r');
    
            xlim([-2 2])
            ylim([-2 2])
            zlim([-2 2])

            hold off
    
            pause(0.1)
        end
        
        figure(2)
        a = -max(max(out.d.Data))/2;
        b = max(max(out.d.Data))/2;
        c = -max(max(out.d.Data))/2;
        d = max(max(out.d.Data))/2;
        e = 0;
        f = max(max(out.d.Data));

        plot3(out.d.Data(:,1), out.d.Data(:,3), out.d.Data(:,2))

        xlim([a b])
        ylim([c d])
        zlim([e f])
        
    else
        disp('plot = false')
    end
    
else
    disp('Error: pry vectors do not exhibit uniform time')
end



function vectarrow(p0,p1,col)
%Arrowline 3-D vector plot.
%   vectarrow(p0,p1) plots a line vector with arrow pointing from point p0
%   to point p1. The function can plot both 2D and 3D vector with arrow
%   depending on the dimension of the input
%
%   Example:
%       3D vector
%       p0 = [1 2 3];   % Coordinate of the first point p0
%       p1 = [4 5 6];   % Coordinate of the second point p1
%       col = 'b'       % Color
%       vectarrow(p0,p1)
%
%       2D vector
%       p0 = [1 2];     % Coordinate of the first point p0
%       p1 = [4 5];     % Coordinate of the second point p1
%       vectarrow(p0,p1)
%
%   See also Vectline
%   Rentian Xiong 4-18-05
%   $Revision: 1.0

  if max(size(p0))==3
      if max(size(p1))==3
          x0 = p0(1);
          y0 = p0(2);
          z0 = p0(3);
          x1 = p1(1);
          y1 = p1(2);
          z1 = p1(3);
          plot3([x0;x1],[y0;y1],[z0;z1],col);   % Draw a line between p0 and p1
          
          p = p1-p0;
          alpha = 0.2;  % Size of arrow head relative to the length of the vector
          beta = 0.1;  % Width of the base of the arrow head relative to the length
          
          hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
          hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
          hw = [z1-alpha*p(3);z1;z1-alpha*p(3)];
          
          hold on
          plot3(hu(:),hv(:),hw(:),'b')  % Plot arrow head
          grid on
          xlabel('x')
          ylabel('y')
          zlabel('z')
          hold off
      else
          error('p0 and p1 must have the same dimension')
      end
  elseif max(size(p0))==2
      if max(size(p1))==2
          x0 = p0(1);
          y0 = p0(2);
          x1 = p1(1);
          y1 = p1(2);
          plot([x0;x1],[y0;y1],col);   % Draw a line between p0 and p1
          
          p = p1-p0;
          alpha = 0.1;  % Size of arrow head relative to the length of the vector
          beta = 0.1;  % Width of the base of the arrow head relative to the length
          
          hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
          hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
          
          hold on
          plot(hu(:),hv(:),col,'b')  % Plot arrow head
          grid on
          xlabel('x')
          ylabel('y')
          hold off
      else
          error('p0 and p1 must have the same dimension')
      end
  else
      error('this function only accepts 2D or 3D vector')
  end
  
end