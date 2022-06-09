function [it, flag, X, err_m] = DLS_multirow(X, toll_err, max_it, alpha, A, W, b,N, n,acc, ptr)

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





while (test==0 && it<max_it && flag==0)
    
    it=it+1;
    if mod(it,500)==0
        it
        err_m
    end
    X_aux=X;
    
    % Update
    for i=1:N
        for j = ptr(i):ptr(i+1)-1
            X(j,i)=(1-alpha)*X(j,i)+alpha*b(j)/A(j,j)-alpha*A(j,[1:j-1,j+1:end])*X([1:j-1,j+1:end],i)/A(j,j);
        end

    end
    
    % Consensus
    
    X_next = zeros(n,N);
    
    for i = 1:N
        X_next(:,i) = sum((repmat(W(i,:),[n,1]).*X),2);
    end
    
    
    % Termination Test
    for i=1:N
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