%% Load nominal robot dynamics parameters

function par = load_wam_links(link_num)

% definitions
ZSFE  =  0.346;              %!< z height of SAA axis above ground
ZHR  =  0.505;              %!< length of upper arm until 4.5cm before elbow link
YEB  =  0.045;              %!< elbow y offset
ZEB  =  0.045;              %!< elbow z offset
YWR  = -0.045;              %!< elbow y offset (back to forewarm)
ZWR  =  0.045;              %!< elbow z offset (back to forearm)
ZWFE  =  0.255;              %!< forearm length (minus 4.5cm)

%% Values taken from Barrett WAM booklet:
% shared/barrett/math/WAM_InertialSpecifications_AC-02.pdf
% link 0 is the base

% link 0 is the base
link0.m = 0.0;
link0.mcm(1) = 0.0;
link0.mcm(2) = 0.0;
link0.mcm(3) = 0.0;
link0.inertia(1,1) = 0.0; 
link0.inertia(1,2) = 0.0; 
link0.inertia(1,3) = 0.0; 
link0.inertia(2,2) = 0.0;   
link0.inertia(2,3) = 0.0;  
link0.inertia(3,3) = 0.0;  

switch link_num 
    case 1
    % SFE joint
    links(1).m = 0.00000; 
    links(1).mcm(1) = -0.00000;
    links(1).mcm(2) = 0.00000;
    links(1).mcm(3) = -0.00000;  
    links(1).inertia(1,1) = 0.00000; 
    links(1).inertia(1,2) = -0.00000; 
    links(1).inertia(1,3) = -0.00000; 
    links(1).inertia(2,2) = 0.00000;   
    links(1).inertia(2,3) = 0.00000;  
    links(1).inertia(3,3) = 0.00000;  
    % SAA joint
    links(2).m = 0.00000;  
    links(2).mcm(1) = 0.43430; 
    links(2).mcm(2) = -0.00495;
    links(2).mcm(3) = -0.00000;
    links(2).inertia(1,1) = 0.24768;  
    links(2).inertia(1,2) = 0.00364;
    links(2).inertia(1,3) = 0.17270;  
    links(2).inertia(2,2) = 0.53601; 
    links(2).inertia(2,3) = -0.02929; 
    links(2).inertia(3,3) =  0.12406;  
    % HR joint    
    links(3).m = 3.53923; 
    links(3).mcm(1) = -0.00889;
    links(3).mcm(2) = -0.02148;
    links(3).mcm(3) = -1.70741;
    links(3).inertia(1,1) = 0.82985;
    links(3).inertia(1,2) = -0.01520;
    links(3).inertia(1,3) = -0.00612;
    links(3).inertia(2,2) = 0.86182;
    links(3).inertia(2,3) = -0.00575;
    links(3).inertia(3,3) = 0.00071;
    % EB joint (elbow)
    links(4).m = 1.03409;
    links(4).mcm(1) = 0.14089;
    links(4).mcm(2) = -0.05914;
    links(4).mcm(3) = -0.00270;
    links(4).inertia(1,1) = 0.01276;
    links(4).inertia(1,2) = 0.00340;
    links(4).inertia(1,3) = -0.00229;
    links(4).inertia(2,2) = 0.02157;
    links(4).inertia(2,3) = 0.00032;
    links(4).inertia(3,3) = 0.03718;
    % WR joint (wrist 1)
    links(5).m = 2.28843;
    links(5).mcm(1) = 0.00709;
    links(5).mcm(2) = 0.00194; 
    links(5).mcm(3) = 0.22347; 
    links(5).inertia(1,1) = 0.02182;
    links(5).inertia(1,2) = -0.00001;
    links(5).inertia(1,3) = -0.00069;
    links(5).inertia(2,2) = 0.02184;
    links(5).inertia(2,3) = -0.00019;
    links(5).inertia(3,3) = 0.00002;
    % WFE joint (wrist 2)
    links(6).m = 0.25655; 
    links(6).mcm(1) = 0.03462; 
    links(6).mcm(2) = -0.00415;  
    links(6).mcm(3) = 0.00121; 
    links(6).inertia(1,1) = 0.00167; 
    links(6).inertia(1,2) = 0.00079; 
    links(6).inertia(1,3) = -0.00009; 
    links(6).inertia(2,2) = 0.00483; 
    links(6).inertia(2,3) = -0.00084;
    links(6).inertia(3,3) = 0.01101;
    % WAA joint (wrist 3)
    links(7).m = 0.63285; 
    links(7).mcm(1) = -0.00157;  
    links(7).mcm(2) = 0.00019;  
    links(7).mcm(3) = 0.08286;  
    links(7).inertia(1,1) = 0.01129; 
    links(7).inertia(1,2) = -0.00006;
    links(7).inertia(1,3) = -0.00006;  
    links(7).inertia(2,2) = 0.01086;
    links(7).inertia(2,3) = 0.00001;
    links(7).inertia(3,3) = 0.00016;
    case 2                     
    % SFE joint
    links(1).m = 10.768; 
    links(1).mcm(1) = links(1).m * -0.0044342;   
    links(1).mcm(2) = links(1).m * 0.12189;    
    links(1).mcm(3) = links(1).m * -0.00066489;
    links(1).inertia(1,1) = 0.29486;    
    links(1).inertia(1,2) = -0.0079502;  
    links(1).inertia(1,3) = -0.00009311;  
    links(1).inertia(2,2) = 0.11350;    
    links(1).inertia(2,3) = -0.00018711;  
    links(1).inertia(3,3) = 0.25065;
    % SAA joint
    links(2).m = 3.8749;
    links(2).mcm(1) = links(2).m * -0.0023698;   
    links(2).mcm(2) = links(2).m * 0.031056;    
    links(2).mcm(3) = links(2).m * 0.015421;
    links(2).inertia(1,1) = 0.026068;  
    links(2).inertia(1,2) = -0.00001346; 
    links(2).inertia(1,3) = -0.00011701;  
    links(2).inertia(2,2) = 0.014722;    
    links(2).inertia(2,3) = 0.00003659;  
    links(2).inertia(3,3) =  0.019348;         
    % HR joint    
    links(3).m = 1.8023;
    links(3).mcm(1) = links(3).m * -0.038259;    
    links(3).mcm(2) = links(3).m * 0.20751;     
    links(3).mcm(3) = links(3).m * 0.00003309;
    links(3).inertia(1,1) = 0.13672;    
    links(3).inertia(1,2) = -0.016804;    
    links(3).inertia(1,3) = 0.00000510;  
    links(3).inertia(2,2) = 0.0058835; 
    links(3).inertia(2,3) = -0.00000530;  
    links(3).inertia(3,3) = 0.13951;
    % EB joint (elbow)
    links(4).m = 2.4007;                        
    links(4).mcm(1) = links(4).m * 0.0049851;
    links(4).mcm(2) = links(4).m * -0.00022942;
    links(4).mcm(3) = links(4).m * 0.13272; 
    links(4).inertia(1,1) = 0.057193;
    links(4).inertia(1,2) = 0.00001467;
    links(4).inertia(1,3) = 0.00008193;
    links(4).inertia(2,2) = 0.057165;
    links(4).inertia(2,3) = -0.00009417;
    links(4).inertia(3,3) = 0.0030044;
    % WR joint (wrist 1)
    links(5).m = 0.12376;                    
    links(5).mcm(1) = links(5).m * 0.00008921;
    links(5).mcm(2) = links(5).m * 0.0051122; 
    links(5).mcm(3) = links(5).m * 0.0043582; 
    links(5).inertia(1,1) = 0.00005587;
    links(5).inertia(1,2) = 0.00000026;
    links(5).inertia(1,3) = 0.00000000;
    links(5).inertia(2,2) = 0.00007817;
    links(5).inertia(2,3) = -0.00000083;
    links(5).inertia(3,3) = 0.00006594;
    % WFE joint (wrist 2)
    links(6).m = 0.41797;                   
    links(6).mcm(1) = links(6).m * 0.00012262; 
    links(6).mcm(2) = links(6).m * -0.017032;  
    links(6).mcm(3) = links(6).m * 0.024683; 
    links(6).inertia(1,1) = 0.00093106; 
    links(6).inertia(1,2) = 0.00000148; 
    links(6).inertia(1,3) = -0.00000201; 
    links(6).inertia(2,2) = 0.00049833; 
    links(6).inertia(2,3) = -0.00022162;
    links(6).inertia(3,3) = 0.00057483;
    % WAA joint (wrist 3)     
    links(7).m = 0.068648; 
    links(7).mcm(1) = links(7).m * -0.00007974;  
    links(7).mcm(2) = links(7).m * 0.00016313;  
    links(7).mcm(3) = links(7).m * -0.0032355;  
    links(7).inertia(1,1) = 0.00003845; 
    links(7).inertia(1,2) = -0.00000019;
    links(7).inertia(1,3) = 0.00000002;  
    links(7).inertia(2,2) = 0.00003878; 
    links(7).inertia(2,3) = -0.00000004;  
    links(7).inertia(3,3) = 0.00007408;
    case 3
    % SFE joint
    links(1).m = 8.3936; 
    links(1).mcm(1) = 0.002942796;   
    links(1).mcm(2) = 0.00528;    
    links(1).mcm(3) = 0.44755;
    links(1).inertia(1,1) = 0.11902;    
    links(1).inertia(1,2) = -0.00010;  
    links(1).inertia(1,3) = -0.00040;  
    links(1).inertia(2,2) = 0.08315;    
    links(1).inertia(2,3) = 0.00068;  
    links(1).inertia(3,3) = 0.29304;
    % SAA joint
    links(2).m = 4.8487;
    links(2).mcm(1) = -0.00108126;   
    links(2).mcm(2) = -0.06485;    
    links(2).mcm(3) = -0.10372;
    links(2).inertia(1,1) = 0.03241;  
    links(2).inertia(1,2) = 0.00011; 
    links(2).inertia(1,3) = -0.00006;  
    links(2).inertia(2,2) = 0.02503;    
    links(2).inertia(2,3) = -0.00274;  
    links(2).inertia(3,3) =  0.20365;         
    % HR joint    
    links(3).m = 1.7251;
    links(3).mcm(1) = 0.01077;    
    links(3).mcm(2) = 0.00004;     
    links(3).mcm(3) = 0.57289;
    links(3).inertia(1,1) = 0.24691;    
    links(3).inertia(1,2) = 0.00001;    
    links(3).inertia(1,3) = -0.00125;  
    links(3).inertia(2,2) = 0.24713; 
    links(3).inertia(2,3) = 0.00000;  
    links(3).inertia(3,3) = 0.07022;
    % EB joint (elbow)
    links(4).m = 2.17266212;                        
    links(4).mcm(1) =  -0.08575;
    links(4).mcm(2) = -0.25915;
    links(4).mcm(3) = 0.00015; 
    links(4).inertia(1,1) = 0.04158;
    links(4).inertia(1,2) = -0.00887;
    links(4).inertia(1,3) = 0.00005;
    links(4).inertia(2,2) = 0.00620;
    links(4).inertia(2,3) = 0.00013;
    links(4).inertia(3,3) = 0.07888;
    % WR joint (wrist 1)
    links(5).m = 0.35655692;                    
    links(5).mcm(1) = 0.00002;
    links(5).mcm(2) = 0.00053; 
    links(5).mcm(3) = 0.09667; 
    links(5).inertia(1,1) = 0.02658;
    links(5).inertia(1,2) = 0.00000;
    links(5).inertia(1,3) = 0.00000;
    links(5).inertia(2,2) = 0.02659;
    links(5).inertia(2,3) = -0.00012;
    links(5).inertia(3,3) = 0.00352;
    % WFE joint (wrist 2)
    links(6).m = 0.40915886;                   
    links(6).mcm(1) = -0.00002; 
    links(6).mcm(2) = -0.00990;  
    links(6).mcm(3) = -0.00690; 
    links(6).inertia(1,1) = 0.00090; 
    links(6).inertia(1,2) = 0.00000; 
    links(6).inertia(1,3) = 0.00000; 
    links(6).inertia(2,2) = 0.00057; 
    links(6).inertia(2,3) = -0.00012;
    links(6).inertia(3,3) = 0.00380;
    % WAA joint (wrist 3)     
    links(7).m = 0.07548270; 
    links(7).mcm(1) = 0.00001;  
    links(7).mcm(2) = 0.00000;  
    links(7).mcm(3) = 0.0042;  
    links(7).inertia(1,1) = 0.00028; 
    links(7).inertia(1,2) = 0.00000;
    links(7).inertia(1,3) = 0.00000;  
    links(7).inertia(2,2) = 0.00028; 
    links(7).inertia(2,3) = 0.00000;  
    links(7).inertia(3,3) = 0.00054;
    otherwise
    error('Link set number not supported!');
end

%% Organize values
% Base position and orientation values are taken from ParameterPool.cf

% make sure inertia matrices are symmetric
for i = 1:7
    for j = 1:3
        for k = j:3
            links(i).inertia(k,j) = links(i).inertia(j,k);
        end
    end
end

% Set default end effector parameters
eff(1).m = 0.0;
eff(1).mcm(1) = 0.0;
eff(1).mcm(2) = 0.0;
eff(1).mcm(3) = 0.0;
eff(1).x(1)  = 0.0;
eff(1).x(2)  = 0.0;
eff(1).x(3)  = 0.30;
eff(1).a(1)  = 0.0;
eff(1).a(2)  = 0.0;
eff(1).a(3)  = 0.0;

% External forces
for j = 1:3
    % I guess this is the external force to the base
    uex0.f(j) = 0.0;
    uex0.t(j) = 0.0;
    for i = 1:7
        uex(i).f(j) = 0.0;
        uex(i).t(j) = 0.0;
    end
end

% base cartesian position and orientation (quaternion)
basec.x  = [0.0,0.0,0.0];
basec.xd = [0.0,0.0,0.0];
basec.xdd = [0.0,0.0,0.0];
baseo.q = [0.0,1.0,0.0,0.0];
baseo.qd = [0.0,0.0,0.0,0.0];
baseo.qdd = [0.0,0.0,0.0,0.0];
baseo.ad = [0.0,0.0,0.0];
baseo.add = [0.0,0.0,0.0];

% Send output
par.basec = basec;
par.baseo = baseo;
par.eff = eff;
par.links = links;
par.link0 = link0;
par.uex = uex;
par.uex0 = uex0;