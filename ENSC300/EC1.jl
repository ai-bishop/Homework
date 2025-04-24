## Extra Credit One:
# Determine the rate of return (internal rate of return) via numerical methods

# initialize values
P = -130000 # initial cost
A = -49000 # annual cost
R = 78000
G = zeros(8,1) # initialize revenue gradient matrix
n = 8
for i = 2:8
    G[i] = G[i - 1] + 1000
end
SV = 23000 # salvage value

# ROR? want i_star
i_star = BreakEven(P, A, R, G, SV, n)




function BreakEven(P, A, R, G, SV, n)

    A_G = A_over_G(G, i, n) # turn gradient of revenue to annuinity
    P_rev = P_over_A(R, i, n) # turn revenue annuinity into present worth
    P_g = P_over_A(A_G, i_star, n) # turn gradient annuinity into present value
    P_A = P_over_A(A, i, n) # turn operating cost into present worth
    P_SV = P_over_F(SV, i, n)

    total = P + P_rev + P_g + P_A + P_SV    

    i = solve(total)
    
    return i
end

function P_over_A(A, i, n)

    P = A * ((1 + i)^n - 1) / ((1 + i) ^ n - 1)

    return P
end

function A_over_G(G, i, n)

    A = G * (1/i - n/((1 + i)^n - 1))

    return A
end

function P_over_F(F, i, n)

    P = F * (1 + i) ^ (-n)
    return P
end
