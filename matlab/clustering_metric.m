function [cluster_cells, g_cc] = clustering_metric(baserun_1)



% THESE HAVE BEEN HARD CODED FOR NOW BUT SHOULD BE CHANGED IF THE DOMAIN  CHANGES
MCDS.mesh.X(1,1) = 10;
MCDS.mesh.X(1,2) = 30;
MCDS.mesh.X(1,end) = 650;
MCDS.mesh.Y(1,1) = 10;
MCDS.mesh.Y(2,1) = 30;
MCDS.mesh.Y(end,1) = 650;

dr = 10; % width of annulus surrounding cells

dom_bounds = [MCDS.mesh.X(1,1)-(MCDS.mesh.X(1,2)-MCDS.mesh.X(1,1))/2, MCDS.mesh.X(1,end)+(MCDS.mesh.X(1,2)-MCDS.mesh.X(1,1))/2; 
              MCDS.mesh.Y(1,1)-(MCDS.mesh.Y(2,1)-MCDS.mesh.Y(1,1))/2, MCDS.mesh.Y(end,1)+(MCDS.mesh.Y(2,1)-MCDS.mesh.Y(1,1))/2]; % domain bounds [x1, x2; y1, y2]
dom_area = (dom_bounds(1,2)-dom_bounds(1,1))*(dom_bounds(2,2)-dom_bounds(2,1)); % whole domain area

R_max = floor(min(dom_bounds(1,2)-dom_bounds(1,1),dom_bounds(2,2)-dom_bounds(2,1))/2); % maximum radius to check - no larger than half of min(dom width,dom height)
num_r = 200; % number of radii to check
r_vals = linspace(0,R_max-dr,num_r); % radii to check

for tindex = 1:61 %this should be the total number of tcounts - leaving as 40 for now. 
    agent_pos = baserun_1.Tcellagent_position_matrix{tindex}(:,1:2); % agent positions (x,y)
    g_cc = zeros(num_r,1); % pair correlation function for T cells to T cells
    num_cells = length(agent_pos);
    for r_n = 1:num_r % for each radius
        r = r_vals(r_n); % get radius
    
        for i = 1:num_cells % for each agent i
    
            % count number of agents in annulus about agent i
            curr_pos = agent_pos(i,:); % agent i's position
            agent_dists = vecnorm(agent_pos-curr_pos,2,2); % distance from current agent to all other agents
            cc_count = 0; % number of cells around cell i
            for j = find(agent_dists > r & agent_dists <= r+dr)' % for all agents j within the annulus surrounding agent i
                cc_count = cc_count + 1; % increment number of cells surrounding cell i
            end
    
            % compute area of annulus about agent i - make sure to account for any area outside domain (assuming it cannot be outside both the left and right at the same time)
            dx = min(abs(curr_pos(1)-dom_bounds(1,1)),abs(curr_pos(1)-dom_bounds(1,2))); % closest distance to left or right edge of domain
            dy = min(abs(curr_pos(2)-dom_bounds(2,1)),abs(curr_pos(2)-dom_bounds(2,2))); % closest distance to top or bottom edge of domain
            innerA = pi*r^2; % area of inner circle (hole of annulus)
            outerA = pi*(r+dr)^2; % area of outer circle (to the outer edge of annulus)
            if dx < r+dr % if the outer circle is poking out of the left or right sides of the domain
                outerA = outerA - ( (r+dr)^2*acos(dx/(r+dr))-dx*sqrt((r+dr)^2-dx^2) ); % subtract the area outside the domain to the left or right from the total area of the outer circle
                if dx < r % if the inner circle is also poking out
                    innerA = innerA - ( r^2*acos(dx/r)-dx*sqrt(r^2-dx^2) ); % subtract the area outside the domain
                end
            end
            if dy < r+dr % if the outer circle is poking out of the top or bottom sides of the domain
                outerA = outerA - ( (r+dr)^2*acos(dy/(r+dr))-dy*sqrt((r+dr)^2-dy^2) ); % subtract the area outside the domain to the top or bottom from the total area of the outer circle
                if dy < r % if the inner circle is also poking out
                    innerA = innerA - ( r^2*acos(dy/r)-dy*sqrt(r^2-dy^2) ); % subtract the area outside the domain
                end
                if dx^2+dy^2 < (r+dr)^2 % if any part of the outer circle is both outside to the left/right and top/bottom (this will only be the case if both dx<r+dr and dy<r+dr)
                    outerA = outerA + (r+dr)^2/2*(acos(dx/(r+dr))+acos(dy/(r+dr))-pi/2) - 1/2*(dx*sqrt((r+dr)^2-dx^2)+dy*sqrt((r+dr)^2-dy^2)) + dx*dy; % add back the corner of area that was subtracted twice from left/right and top/bottom
                    if dx^2+dy^2 < r^2 % if the inner circle is also outside the left/right and top/bottom
                        innerA = innerA + r^2/2*(acos(dx/r)+acos(dy/r)-pi/2) - 1/2*(dx*sqrt(r^2-dx^2)+dy*sqrt(r^2-dy^2)) + dx*dy; % add back the piece that was subtracted twice
                    end
                end
            end
    
            % computer density of agents in annulus
            annulusA = outerA - innerA; % area of annulus is then outer circle minus inner circle
            g_cc(r_n) = g_cc(r_n) + cc_count/annulusA; % add density of cells around cell i to the total
        end
    end
    g_cc = g_cc*dom_area/num_cells^2; % average density of T cells within radius r of T cells, divided by mean-field density of T cells across whole domain
    g_cc_store{tindex} = g_cc;

    cluster_cells(tindex) = trapz(r_vals, g_cc)/(R_max-dr); % average value of pair correlation function to get measure of clustering

end
end