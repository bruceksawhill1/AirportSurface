%  Script to compute constrained permutations of k integers
twoperms = cell(10,1); % create cell to store permutations
for k = 1: 10
    
legals = [];
testperms = perms(1:k);  % Start with all possible permutations
pl = length(testperms);
pw = length(testperms(1,:));

for i = 1: length(testperms); % Consraint implementation loop
    
    productindex = 1;
    for j = 1: length(testperms(1,:))
     % this if statement encodes constraint (no more than n positions
     % out of place, encoded as >= n + 1
        if abs(testperms(i,j)-j) >= 3, % test out-of-place-ness
          productindex = 0;
        end
        
    end
    
    if productindex == 1
         legals = [legals; testperms(i,:)];
    end
    
end
twoperms{k} = legals;
end

legals
length(legals)