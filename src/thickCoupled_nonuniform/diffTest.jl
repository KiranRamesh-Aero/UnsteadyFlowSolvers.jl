


dqdthe = diff(qu)./diff(surf.theta)

dthedx = (2/surf.c)*sin.(surf.theta)

dqdx = dqdthe.*dthedx[1:end-1]
