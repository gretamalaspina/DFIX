function [it, flag, X, err_m] = DLS(X, toll_err, max_it, alpha, A, W, b, n,acc)

% X is a nxn matrix, X(:,i) is the iterate in the i-th node

% toll_err : tolerance for the termination test
% max_it : maximum number of iterations
% alpha : stepsize
% A : coeficient matrix of the system
% b : RHS of the system
% W : weights matrix associated with the graph
% n : dimension

% flag = 1 if the termination condition is satisfied within max_it iterations
%       -1 otherwise
% X is a nxn matrix, X(:,i) is the iterate in t he i-th node
flag=0;
it = 0;
test = 0;

if acc==1
    EV=eig(W);
    sigma=EV(end-1);
    g=(1-sqrt(1-sigma^2))/(1+sqrt(1-sigma^2));
else
    g=0;
end


while (test==0 && it<max_it && flag==0)
    
    it=it+1;
    if mod(it,500)==0
        it
        err_m
    end
    X_aux=X;
    
    % Update
    for i=1:n
        X(i,i)=(1-alpha)*X(i,i)+alpha*b(i)/A(i,i)-alpha*A(i,[1:i-1,i+1:end])*X([1:i-1,i+1:end],i)/A(i,i);
    end
    
    % Consensus
    
    X_next = zeros(n,n);
    for i = 1:n
        X_next(:,i) = (1+g)*sum((repmat(W(i,:),[n,1]).*X),2)-g*X_aux(:,i);
    end
    
    
    % Termination Test
    for i=1:n
        err(i) = norm(A*X_next(:,i)-b);
    end
    err_m = max(err);
    if err_m<toll_err
        test = 1;
        flag = 1;
    end
    
    X = X_next;
    
    
end



if it==max_it
    flag = -1;
    
end

end