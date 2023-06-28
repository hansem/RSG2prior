function angle=estAngle(v,w)

% estimate angle between two vectors
% v and w should in the same size

angle=acosd(dot(normalize(v(:),1),normalize(w(:),1),1));