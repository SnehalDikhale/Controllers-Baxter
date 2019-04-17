function [k,p] = dynamics(cm,m)

 syms t1(t) t2(t) t3(t) t4(t) t5(t) t6(t) t7(t); 
 syms q1 q2 q3 q4 q5 q6 q7 g ;

% Substituting symbolic theta variables (as a function of time) into 
% center of mass expressions.
cm = subs(cm,[q1,q2,q3,q4,q5,q6,q7],[t1(t),t2(t),t3(t),t4(t),t5(t),...
        t6(t),t7(t)]);
      v = diff(cm, t);
     k = (1/2)*m*(transpose(v)*v);
     p = m * g * (transpose(cm) * [0; 0; 1]);