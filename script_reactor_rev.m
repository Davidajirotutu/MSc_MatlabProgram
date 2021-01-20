%% Solving Nuclear Diffusion in 3D 
% Import 3D Reactor Geometry
tic
model = createpde(2);
importGeometry(model,'model_IAEA_340h.STL');
pdegplot(model,'FaceLabels','on');
%%  and inputs all Data
rectarray=[0 0 0 0 0 0 4 0 0;
           0 0 0 0 0 1 4 4 0;
           0 0 0 0 3 1 1 4 0;
           0 0 0 2 2 2 1 4 4;
           0 0 2 2 2 2 1 1 4;
           0 2 2 2 2 2 2 1 4;
           3 2 2 2 3 2 2 1 4];
nonfuel= 4;  % Non fuel lattice       
max_row = size(rectarray,1);
max_column = size(rectarray,2);
r=[0;0]; % r is for dirichlet bc
g=[0;0]; % g is all zeros for both reflective and vacuum bc
q_ref= zeros(2); % q is all zeros for reflective bc
q_ref=q_ref(:); % q_ref turned to a column vector
q_vac=(0.46922)*eye(2); % q_vac is all 0.4692 diagonal matrix for reflective bc
q_vac=q_vac(:); % q_vac turned to a column vector
coreshape=45; % 1/8 model
lengtharray_x= [0.5 1 1 1 1 1 1 1 1 1]*20;
lengtharray_y= [0.5 1 1 1 1 1 1 1]*20;
lengtharray_z = zeros(1,10)*20;
length_z= 380;
llcoordinate=[0 0 0]; 
%% Calculate the absorption with the axial bulking
axial_bucking = 0.8e-4;
D= [1.5 0.4;1.5 0.4;1.5 0.4;2.0 0.3;2.0 0.3];
ab=[0.01 0.08;0.01 0.085;0.01 0.13;0.0 0.01;0.0 0.055]; 
%ab = axial_bucking*D + absorption;

% define cross sections
xs = struct('D',{},'a',{},'s',{},'vf',{},'f',{},'X',{});
xs(1).D=D(1,:);   xs(2).D=D(2,:);     
xs(3).D=D(3,:);    xs(4).D=D(4,:);    
xs(5).D=D(5,:); 

xs(1).a=ab(1,:); xs(2).a=ab(2,:);
xs(3).a=ab(3,:); xs(4).a=ab(4,:); 
xs(5).a=ab(5,:);

xs(1).f=[0.0 0.0];  xs(2).f=[0.0 0.0];  
xs(3).f=[0.0 0.0];  xs(4).f=[0.0 0.0]; 
xs(5).f=[0.0 0.0];

xs(1).vf=[0.0 0.135];  xs(2).vf=[0.0 0.135];  
xs(3).vf=[0.0 0.135];  xs(4).vf=[0.0 0.0];
xs(5).vf=[0.0 0.0];

xs(1).s=[0 0.02; ... %scattering 1->1, 1->2
         0 0];       %scattering 2->1, 2->2. g->g is not used in diffusion theory. Enter any value like zero.
xs(2).s=[0 0.02; 0 0]; xs(3).s=[0 0.02; 0 0]; 
xs(4).s=[0 0.04; 0 0]; xs(5).s=[0 0.04; 0 0]; 

xs(1).X=[1 0]; xs(2).X=[1 0]; xs(3).X=[1 0]; xs(4).X=[1 0]; xs(5).X=[1 0];
%% Apply Boundary Condition
applyBoundaryCondition(model,'neumann','Face',3:7,'g',g,'q',q_vac);
%% Specify Coefficient & Apply Boundary condition & Solve 
[xcoordinates,ycoordinates] = rectarraymul_llcoordinates1(rectarray,lengtharray_x,lengtharray_y,llcoordinate);
 m = 0;
 f = 0;
d= @(location, ~) dfission(xs,xcoordinates,ycoordinates,rectarray,location);
a= @(location, ~) aabsorption(xs,xcoordinates,ycoordinates,rectarray,location);
c= @(location, ~) cdiffusion(xs,xcoordinates,ycoordinates,rectarray,location);
specifyCoefficients(model,'m',m,'d',d,'c',c,'a',a,'f',f);
generateMesh(model,'Hmax',7.0)
evr=[0 1.5];
result = solvepdeeig(model,evr);
plotpdeeig(model,result)
%% Calculate the fission reaction rate
totalpower = 22.125; %total_fuel(rectarray,nonfuel);
power_array= power_rect_3DD(result,rectarray,llcoordinate,lengtharray_x,lengtharray_y,xs,totalpower);
toc