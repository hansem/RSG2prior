function w=randVect(v,a,n)

% w=randVect(v,a,n)
% generate n random unit vectors w whose angle with vector v is a
% input: v [nD x 1],a,n
% output: w [nD x n]
% algorithm
% from CrossValidated: "generating vectors under constraints on 1 and 2 norm"
% 1) random vector x on (n-1) dimensional unit sphere> [0 x]
% 2) orthgonal matrix U whose 1st row is v (use GramSchmidt)
% 3) r=sqrt(1-(cosd(a)^2)/(v'*v)), x2=r*U'*[0;x(:)] (norm(x2,2)=r, (v(:)')*x2=0)
% 4) w=x2+cosd(a)/(v'*v)*v; (check norm(w,2) and v(:)'*w)

% depends on GramSchmidt.m

% old
% 1) find hyperplane whose |angle with v| < a+eps: dot(v,w)=w1*v1+...+wn*vn=cos(a)
% 2) sample s.t. unity constraints: sum(w.^2)=1

%% init
nD=length(v);
v=normalize(v(:),1);
i=1;

err=eps; % (-3); % acceptable error in angle; eps~10^(-16)

try 
    w=nan(nD,n);
catch
    disp('too many samples requested: out of memory');
end

%% main
while i<=n     
    % 1) random vector x on (n-1) dimensional unit sphere> [0 x]
    x=[0; normalize(randn(nD-1,1),1)]; x=x(:)';
    
    % 2) orthgonal matrix U whose 1st row is v (use GramSchmidt)
    U=[v rand(nD,nD-1)];
    U=GramSchmidt(U); % orth use SVD
    U=U';
    
    % 3) r=sqrt(1-(cosd(a)^2)/(v'*v)), x2=r*U'*[0;x(:)] (norm(x2,2)=r, (v(:)')*x2=0)
    r=sqrt(1-(cosd(a)^2)/(v'*v));
    x2=r*U'*x(:);
%     disp(r);
%     disp(norm(x2,2));
%     disp((v(:)')*x2);
    
    % 4) w=x2+cosd(a)/(v'*v)*v; (check norm(w,2) and v(:)'*w)
    wTmp=x2+cosd(a)/(v'*v)*v;
%     disp(norm(wTmp,2));
%     disp(v(:)'*wTmp);
    if abs(norm(wTmp,2)-1)<err & abs(v(:)'*wTmp-cosd(a))<err
        w(:,i)=wTmp;
        i=i+1;
    end
    
    % old
%     % 1) find hyperplane whose |angle with v| < a+eps: dot(v,w)=w1*v1+...+wn*vn=cos(a)
%     % idea: sample uniform random integer>sort>difference>divide by their sum (Kraemer algorithm)
%     x=randi([0 2^52],nD+1,1); % intmax('uint64')
%     y=diff(sort(x));
%     y=y./sum(y);
%     
%     wTmp=y*cosd(a)./v;
% %     disp([num2str(angleVectors(v,wTmp)) ', ' num2str(norm(wTmp))]);
%     
%     % 2) sample s.t. unity constraints: sum(w.^2)=1
%     if abs(norm(wTmp)-1)<err
%         w(:,i)=wTmp;
% %         if rem(i,round(n/10))==0
%             disp([num2str(i) '/' num2str(n) ' done']);
% %         end        
%     end % if abs(norm(wTmp)-1)<err
    
%     % old: 2)>1)
%     rv=normalize(randn(nD,1),1);
%         
%     if abs(acosd(dot(rv,v))-a)<=err
%         w(:,i)=rv;
%         
% %         if rem(i,round(n/10))==0
%             disp([num2str(i) '/' num2str(n) ' done']);
% %         end
%         i=i+1;
%     end % if abs(angleVectors(rv,v)-a)<=err    
end % while i<=n

%% test
function test
v=randn(10,1);
a=10;
n=100;
tic;
w=randVect(v,a,n);
toc;
% check angle histogram
angle=[];
for i=1:size(w,2)
    angle=[angle;angleVectors(v,w(:,i))];
end
figure; hist(angle,100);
% check vectors are different
d=zeros(n,n);
for i=1:(n-1)
    for j=(i+1):n
        d(i,j)=sqrt(sum((w(:,i)-w(:,j)).^2));
    end
end
figure; imagesc(d);