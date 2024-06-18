colorVec = [[0,0.4470,0.7410];        %blueish (good for plotting)
            %[0.8500, 0.3250, 0.0980]; %reddish
            [0.85, 0.33, 0.1];        %orangeish (better than the one above?)
          	[0.9290, 0.6940, 0.1250]; %yellowish (good for plotting though)
          	[0.4940, 0.1840, 0.5560]; %purple
          	[0.4660, 0.6740, 0.1880]; %light green
          	[0.3010, 0.7450, 0.9330]; %light blue
          	[0.6350, 0.0780, 0.1840]; %maroon/brown
            1/255*[0,104,87];         %dark green
            1/255*[0,102,0];          % also dark green
            1/255*[153,0,0];          % also maroon/brown
            1/255*[204,102,0];        %dark orange
            1/255*[102,0,102];        %dark purple
            1/255*[0,0,153];          %dark blue
            1/255*[255,0,127];        %reddish pink
            1/255*[160,160,160];      %gray
            [0,0,0]                   %black
           ];
% indices 1, 2, 5 are good together for plotting. Also 3 and 4. 6 too, if
% necessary. 16 is black.
iiBlack = 16;