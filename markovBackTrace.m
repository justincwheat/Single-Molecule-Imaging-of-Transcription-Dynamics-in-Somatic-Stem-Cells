function B = markovBackTrace(A,v0,t)
%Input: A = the Markov Chain
%       v0 = the initial distribution
%       t = the elapsed "time" (number of cell divisions)
%
%Output: B = a matrix whose ijth coordinate is the expected number of times
%            the trajectory occupied state j given it occupies state i at
%            time t.
%
%
n=length(A);
B=zeros(n);
vt=(A^t)*v0;
for i=1:n
    for j=1:n
        for s=0:(t-1)
            vs=(A^s)*v0;
            Ainter=A^(t-s);
            B(i,j)=B(i,j)+Ainter(i,j)*vs(j);
        end
        B(i,j)=B(i,j)/vt(i);
    end
end

end