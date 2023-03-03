function x_best = differential_evolution(func, F, CR, NP, D, Gmax)
% MATLAB implementation of the Differential Evolution for cost function
% minimization


x = zeros(D, NP);
v = zeros(D, NP);
u = zeros(D, NP);
x2 = zeros(D, NP);

score = zeros(1,NP);

ql = -1; % minimum quaternion element value
qu = 1; % maximum quaternion element value
wl = -2*pi/180; % minimum angular rate
wu = 2*pi/180; % maximum angular rate
ml = -0.05; % minimum dipole
mu = 0.05; % maximum dipole
bl = -4*10^3; %minimum bias
bu = 4*10^3; %maximum bias

% rng('default');

%-----------------------------------------------------------------------------
%% Initialization
%-----------------------------------------------------------------------------

for i=1:NP
    
    %for k=1:D
    randoms = rand([13,1]);
    
    x(:,i) = [ql + randoms(1)*(qu - ql);
              ql + randoms(2)*(qu - ql);
              ql + randoms(3)*(qu - ql);
              ql + randoms(4)*(qu - ql);
              wl + randoms(5)*(wu - wl);
              wl + randoms(6)*(wu - wl);
              wl + randoms(7)*(wu - wl);
              ml + randoms(8)*(mu - ml);
              ml + randoms(9)*(mu - ml);
              ml + randoms(10)*(mu - ml);
              bl + randoms(11)*(bu - bl);
              bl + randoms(12)*(bu - bl);
              bl + randoms(13)*(bu - bl)];
          
    x(1:4,i) = x(1:4,i)/norm(x(1:4,i));
    %end 
end


% COMPUTE COST FOR EACH VECTOR
for i=1:NP
    score(i) = func(x(:,i));
end

% Choose the best score = minimum squared error
[best_value, idx] = min(score);

G = 0;


%-----------------------------------------------------------------------------
%% Mutation, Crossover and Selection
%-----------------------------------------------------------------------------


while(G <= Gmax)
    
    
    for i=1:NP
       
        in_area_w = 0; % A flag defining if the angular velocity of the mutated vector is in acceptable area
        in_area_m = 0;
        in_area_b = 0;
        
        % MUTATION
        
        % Indexes must be different
        while(in_area_w == 0 || in_area_m == 0 || in_area_b == 0)
            r1 = randi([1 NP]);
            while(i == r1)
                r1 = randi([1 NP]);
            end
            
            r2 = randi([1 NP]);
            while(i==r2)||(r2==r1)
                r2 = randi([1 NP]);
            end
            
            r3 = randi([1 NP]);
            while((i==r3)||(r3==r1)||(r3==r2))
                r3 = randi([1 NP]);
            end
            
            %  indexes must be different
    %         while((i==r1)||(i==r2)||(i==r3))
    %             r1 = randi([1 NP]);
    %             r2 = randi([1 NP]);
    %             r3 = randi([1 NP]);
    %         end
            
            % Donor vector
            v(:,i) = x(:, r1) + F*(x(:, r2) - x(:, r3));
            
            % If the angular velocity components are more than wu, the mutation repeats
            if max(abs([v(5,i),v(6,i),v(7,i)]))<=wu 
                in_area_w = 1; % The angular velocity is in acceptable area
            end
            
            if max(abs([v(8,i),v(9,i),v(10,i)]))<=mu
                in_area_m = 1;
            end

            if max(abs([v(11,i),v(12,i),v(13,i)]))<=bu
                in_area_b = 1;
            end
        end
    
        v(1:4,i) = v(1:4,i)/norm(v(1:4,i));
        
        
        
        % CROSSOVER
        jrand = randi([1 D]);
        for j=1:D

            rd = rand();
            %jrand = randi([1 D]);

            if (rd <= CR)||(j == jrand)
                u(j,i) = v(j,i);
            else
                u(j,i) = x(j,i);
            end

        end % FOR

        u(1:4,i) = u(1:4,i)/norm(u(1:4,i));

        
        
        % SELECTION
        u_cost = func(u(:,i));

        if u_cost<=score(i)
            x2(:,i) = u(:,i);
            score(i) = u_cost;
        else
            x2(:,i) = x(:,i);
        end
    
    end % FOR
    
    
    x = x2;
    
    G = G+1;
    
%     for i=1:NP
%         score(i) = func(x(:,i));
%     end
    
    best_old = best_value;
    idx_old = idx;
    
    [best_value, idx] = min(score);
    best_value_plot(G) = best_value;
    
    disp(G);
    disp(best_value);
end % WHILE


%-----------------------------------------------------------------------------
%% Plot of the cost function value by number of generations
%-----------------------------------------------------------------------------


figure
plot([0:1:Gmax],best_value_plot,'LineWidth',2)
xlabel('Number of generation')
ylabel('Cost function value')
grid on

x_best = x(:,idx);

end