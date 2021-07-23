function [positions]=Constrained_Random_Walk(wall_thickness, opening_size, step_size, N_steps, N_cells)
    %{
    Lets code a constrained random walk that acts as a hidden markov
    process for biased diffusion. That is, when it is on the right half of
    the box, the cell is biased or persistent, but when on the left side of the box the
    cell is random.

    Some trajectories may APPEAR to go across the impermeable wall, but
    if you check the length, it is ALWAYS less than step_size. So what
    happens is a reflection first followed by propogation, then a random step. The line that is drawn
    connects the beginning and end, but it doesn't "actually" go through
    the wall, it just doesn't draw the midpoint.

    %}
    %the only question is what happens when the particle hits the edge in the
    %random walk? Does it stop? Does it reflect? We're doing a reflect. 
    Box_Width = 10;
    %pick a random direction 
    theta = rand(N_cells,N_steps)*2*pi;
    wall_thickness = wall_thickness*Box_Width;
    opening_size = opening_size*Box_Width;

    %how do we decide reflective positions?
    %well, anything <0 or >Box_Width is reflective. 
    %also, anything from x = (Box_Width / 2) +- (wall_thickness/2)
    %for all y except y = (Box_Width / 2) +- (opening_size / 2)

    wall_start = (Box_Width/2)-(wall_thickness/2);
    wall_end = (Box_Width/2)+(wall_thickness/2);
    opening_start = (Box_Width/2)-(opening_size/2);
    opening_end = (Box_Width/2)+(opening_size/2);
    %pick a random point to start that is within the box. 
    start_positions = rand(N_cells,2)*Box_Width;

    C = 1;
    while C==1
        %indeces in lower wall
        xbad = (start_positions(:,1)>wall_start).*(start_positions(:,1)<wall_end);
        ybad = start_positions(:,2)<opening_start;
        bad_inds = (xbad.*ybad) > 0;
        %indeces in upper wall
        ybad = start_positions(:,2)>opening_end;
        bad_inds((xbad.*ybad)>0)=true;
        if sum(bad_inds)==0
            C=0;
        else
            start_positions(bad_inds,:) = rand(sum(bad_inds),2);
        end
    end

    positions = zeros(N_cells,2,N_steps);
    positions(:,:,1)=start_positions;

    %equations of boundaries?
    %they form line segments.
    %Segment = {(X1, Y1), (X2, Y2)}
    bottom_left_line = [0,0, wall_start,0];

    left_bottom_wall = [wall_start,0, wall_start,opening_start];
    top_bottom_wall = [wall_start,opening_start, wall_end,opening_start];
    right_bottom_wall = [wall_end,opening_start, wall_end,0];

    bottom_right_line = [wall_end,0, Box_Width,0];
    right_right_line = [Box_Width,0, Box_Width,Box_Width];
    top_right_line = [Box_Width,Box_Width, wall_end,Box_Width];

    right_top_wall = [wall_end,Box_Width, wall_end,opening_end];
    bottom_top_wall = [wall_end,opening_end, wall_start,opening_end];
    left_top_wall = [wall_start,opening_end, wall_start,Box_Width];

    top_left_line = [wall_start,Box_Width, 0, Box_Width];
    left_left_line = [0, Box_Width, 0,0];

    boundaries = {bottom_left_line, left_bottom_wall, top_bottom_wall, right_bottom_wall,...
                  bottom_right_line, right_right_line, top_right_line, right_top_wall, ...
                  bottom_top_wall, left_top_wall, top_left_line, left_left_line};

    parallels = cell(1,length(boundaries)); %runs counterclockwise starting in bottom left. 

    %normals = {[0,1],[-1,0],[0,1],[1,0],[0,1],[-1,0],[0,-1],[1,0],[0,-1],[-1,0],[0,-1],[1,0]};
    %normal vectors. Realistically, we should get these by rotating boundary
    %segment.
    normals = cell(1,length(boundaries));
    %rotation matrix (90 degrees CCW)
    R = [cos(pi/2), -sin(pi/2) ; sin(pi/2), cos(pi/2)];

    for seg = 1:length(parallels)
        %parallels is really just = [deltax, deltay] / norm([deltax, deltay])
        deltax = boundaries{seg}(3) - boundaries{seg}(1);
        deltay = (boundaries{seg}(4) - boundaries{seg}(2));
        parallels{seg} = [deltax,deltay]/norm([deltax,deltay]);

        normals{seg} = (R * parallels{seg}')';
    end


    %figure(5)
    %hypothetically, if you plot all those, you should see your boundary. 
    % for b=1:length(boundaries)
    %     xs = [boundaries{b}(1), boundaries{b}(3)];
    %     ys = [boundaries{b}(2), boundaries{b}(4)];
    %     plot(xs,ys,'-k');
    %     hold on
    % end
    %looks good. 
    %plot the starting positions with red x's.
    % plot(positions(:,1,1),positions(:,2,1),'rx')
    % hold on
    %next, need to simulate the steps.
    %how do you take care of if it intersects the wall?

    wb = waitbar(0,'Progress shown here');
    wb_i=1;
    for t = 2:N_steps
        wb_i=wb_i+1;
        steps = [step_size*cos(theta(:,t)), step_size*sin(theta(:,t))];
        positions(:,:,t)=positions(:,:,t-1)+steps;

        positions_f = positions(:,:,t);
        positions_past = positions(:,:,t-1);
        %plot them.
        %     for n = 1:N_cells
        %         xs = positions(n,1,t-1:t);
        %         ys = positions(n,2,t-1:t);
        %         plot(xs(:),ys(:),'-rx')
        %         hold on
        %         plot(xs(1),ys(1),'gx');
        %         hold on
        %     end
        time_steps = ones(N_cells,1);

        check=1;
        bad_inds = true(N_cells,1);

        while check==1
            %so bad_inds here are a red flag. 
            %for these, determine which boundary segments they intersected.

            %so each makes a line segment with from (xi, yi, xf, yf)
            %(positions(:,1,t-1), positions(:,2,t-1), positions(:,1,t), position(:,2,t))
            %for each line segment, check if they have mutual interval existance. 

            %let S1 = boundary = [x1, y1, x2, y2]
            %let S2 = step = [x3, y3, x4, y4]
            %they have mutual overlap if 
            %min([x1,x2])<max([x3,x4] AND max([x1,x2])>min([x3,x4])
            %AND min([y1,y2])<max([y3,y4]) AND max([y1,y2])>min([y3,y4])

            %if they intersect, you can determine where using
            %yi = m_1*xi + b_1 --> m_1 = slope , b_1 = y-intercept
            %yi = m_2*xi + b_2
            %x_i = (b_1 - b_2) / (m_2 - m_1)
            %y_i = (m_1 * (b_1 - b_2) / (m_2 - m_1)) + b_1

            if sum(bad_inds)<1
                check=0;
            end

            useful_inds = find(bad_inds>0);
            for ind = 1:sum(bad_inds)
                bad_i = useful_inds(ind);

                overlap_bounds = [];
                loc_intercept = {};
                
                x1 = positions_past(bad_i,1);
                y1 = positions_past(bad_i,2);
                x2 = positions_f(bad_i,1);
                y2 = positions_f(bad_i,2);
                
                for seg = 1:length(boundaries)
                    x3 = boundaries{seg}(1);
                    x4 = boundaries{seg}(3);
                    y3 = boundaries{seg}(2);
                    y4 = boundaries{seg}(4);

                    if min([x1,x2])<max([x3,x4]) && max([x1,x2])>min([x3,x4]) && min([y1,y2])<max([y3,y4]) && max([y1,y2])>min([y3,y4])
                        %overlap! Potentially. This isn't overlap in special
                        %cases

                        %where does it overlap?
                        %calculate slope of line velocity and slope of boundary.
                        m_boundary = (y4-y3)/(x4-x3);
                        m_line = (y2-y1)/(x2-x1);

                        if isinf(m_line) && m_boundary~=0
                            %vertical line segment, sloping boundary
                            xi = x2;
                            yi = (m_boundary*(xi - x3)) + y3;
                        elseif m_line==0 && ~isinf(m_boundary)
                            %horizontal line segment, sloping boundary
                            yi = y1;
                            xi = ((yi - y3)/m_boundary) + x3;

                        elseif isinf(m_boundary) && m_line~=0
                            %vertical boundary, sloping line segment
                            xi = x3;
                            yi = (m_line*(xi -x1)) + y1;
                        elseif m_boundary==0 && m_line~=0
                            %horizontal boundary, sloping line segment
                            yi = y3;
                            xi = ((yi-y1)/m_line) + x1;
                        elseif isinf(m_boundary) && m_line==0
                            %vertical boundary, horizontal line segment. 
                            xi = x3;
                            yi = y1;
                        elseif m_boundary==0 && isinf(m_line)
                            %horizontal boundary, vertical line segment
                            xi = x1;
                            yi = y3;
                        else
                            %both boundary and line segment are sloping.
                            %note, there are no parallel situations, since they
                            %wouldn't intercept.
                            %determine y intercepts of equations.
                            b_line = y1-(m_line*x1);
                            b_boundary = y3-(m_boundary*x3);
                            xi = ( ((-m_boundary*x3)+b_boundary) + ((m_line*x1)-b_line) ) / (m_line - m_boundary);
                            yi = (m_line*(xi - x1)) + b_line + y1;
                        end
                        %does xi exist between x1 and x2, and between x3 and x4?
                        %between_line_x = 
                        %does yi exist between y1 and y2, and between y3 and y4?
                        %if xi<x2 && xi>x1 && xi<x4 && xi>x3 %not right, could
                        if xi<=max([x1,x2]) && xi>=min([x1,x2]) && xi<=max([x3,x4]) && xi>=min([x3,x4])
                            %x is good
                            if yi<=max([y1,y2]) && yi>=min([y1,y2]) && yi<=max([y3,y4]) && yi>=min([y3,y4])
                                overlap_bounds(end+1) = seg;
                                loc_intercept{end+1} = [xi,yi];
                            end
                        end
                    end
                end
                %if multiple intercept locations, check to make sure they are the
                %same.
                if length(loc_intercept)>1
                    %determine time of impacts. 
                    Im_t = zeros(length(loc_intercept),1);
                    vx = (x2 - x1)/time_steps(bad_i);
                    vy = (y2 - y1)/time_steps(bad_i);
                    for seg_i = 1:length(loc_intercept)
                        deltax=loc_intercept{seg_i}(1) - x1;
                        deltat = abs(deltax / vx);
                        Im_t(seg_i) = deltat;
                    end
                    %determine first impact
                    [min_val,min_ind] = min(Im_t);
                    if sum(Im_t==min_val)>1
                        %if impacted corner? They would be equal.
                        %should have impacted corner, should only be length 2
                        %determine time of collision = deltat
                        deltat = min_val;
                        %reflection direction is opposite of velocity
                        vf = [-vx,-vy];
                    else
                        %did not impact corner. 
                        %use min_ind 
                        deltat = Im_t(min_ind);
                        %reflection direction needs to be calculated
                        B = -dot(normals{overlap_bounds(min_ind)},[vx,vy])*normals{overlap_bounds(min_ind)};
                        A = dot(parallels{overlap_bounds(min_ind)},[vx,vy])*parallels{overlap_bounds(min_ind)};
                        vf = A+B;

                    end
                    xi = loc_intercept{min_ind}(1);
                    yi = loc_intercept{min_ind}(2);
                    %propogate f for 1-deltat time.
                    %should actually be time_step - deltat
                    t_remain = time_steps(bad_i)-deltat;
                    xf = xi + vf(1)*t_remain;
                    yf = yi + vf(2)*t_remain;
                    %update positions and time
                    positions_f(bad_i,:) = [xf,yf];
                    positions_past(bad_i,:) = [xi,yi];
                    %plot([xi,xf],[yi,yf],'-bx')
                    %hold on
                    time_steps(bad_i)=t_remain; 

                elseif length(loc_intercept)==1
                    %determine the time of collision. 
                    %determine the reflection direction.
                    %propogate f
                    xi = loc_intercept{1}(1);
                    yi = loc_intercept{1}(2);
                    %determine time of collision = deltat
                    deltax=xi-positions_past(bad_i,1);
                    vx = (positions_f(bad_i,1) - positions_past(bad_i,1))/time_steps(bad_i);
                    vy = (positions_f(bad_i,2) - positions_past(bad_i,2))/time_steps(bad_i);
                    deltat = abs(deltax / vx);
                    %reflection direction needs to be calculated
                    B = -dot(normals{overlap_bounds(1)},[vx,vy])*normals{overlap_bounds(1)};
                    A = dot(parallels{overlap_bounds(1)},[vx,vy])*parallels{overlap_bounds(1)};
                    vf = A+B;
                    %propogate f for 1-deltat time.
                    %should actually be time_step - deltat
                    t_remain = time_steps(bad_i)-deltat;
                    xf = xi + vf(1)*t_remain;
                    yf = yi + vf(2)*t_remain;
                    %plot([xi,xf],[yi,yf],'-bx')
                    %hold on
                    %update positions and time
                    positions_f(bad_i,:) = [xf,yf];
                    positions_past(bad_i,:) = [xi,yi];
                    time_steps(bad_i)=t_remain; 
                else
                    %disp("no intersection")
                    %plot([positions_past(bad_i,1),positions_f(bad_i,1)],[positions_past(bad_i,2),positions_f(bad_i,2)],'-bx')
                    %hold on
                    bad_inds(bad_i)=0;

                end
            end

        end %end while loop 

        %replace final positions
        %plot corrections to check
        %     for n = 1:N_cells
        %         xs = [positions(n,1,t-1), positions_f(n,1)];
        %         ys = [positions(n,2,t-1), positions_f(n,2)];
        %         plot(xs(:),ys(:),'-kx')
        %         hold on
        %     end
        positions(:,:,t) = positions_f;
        waitbar(wb_i/N_steps,wb)

        %     clf(5) %clear figure
        %     %replotting
        %     for b=1:length(boundaries)
        %         xs = [boundaries{b}(1), boundaries{b}(3)];
        %         ys = [boundaries{b}(2), boundaries{b}(4)];
        %         plot(xs,ys,'-k');
        %         hold on
        %     end
        %     for n = 1:N_cells
        %         xs = positions(n,1,1:t);%[positions(n,1,t-1), positions_f(n,1)];
        %         ys = positions(n,2,1:t);%[positions(n,2,t-1), positions_f(n,2)];
        %         plot(xs(:),ys(:),'-kx')
        %         hold on
        %         plot(xs(end),ys(end),'gx')
        %         hold on
        %     end

    end
    figure(5000)
    %plot boundaries
    for b=1:length(boundaries)
        xs = [boundaries{b}(1), boundaries{b}(3)];
        ys = [boundaries{b}(2), boundaries{b}(4)];
        plot(xs,ys,'-k');
        hold on
    end
    %plot trajectories
    for n=1:min([10,N_cells])
        xs = positions(n,1,:);
        ys = positions(n,2,:);
        plot(xs(:),ys(:),'-x')
        hold on
    end
    close(wb)
end
%so need to find the intersection point of 
