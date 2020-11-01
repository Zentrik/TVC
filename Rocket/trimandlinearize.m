xKeep = {...
    'q0 q1 q2 q3(1)'
    'q0 q1 q2 q3(2)'
    'q0 q1 q2 q3(3)'
    'q0 q1 q2 q3(4)'
    'ub,vb,wb(1)'
    'ub,vb,wb(2)'
    'p,q,r(1)'
    'p,q,r(2)'
    'p,q,r(3)'};

[~, xElim] = setdiff(G7.StateName, xKeep);
