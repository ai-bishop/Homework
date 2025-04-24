% calc cash flows
cash_flow = [-130000, 29000, 30000, 31000, 32000, 33000, 34000, 35000, 59000];

% create sum of cash flow as func of r, as present values
cf_sum = @(r) sum(cash_flows ./ (1 + r).^(0:length(cash_flows)-1));

% initial guess for i_star
i_init = 0.1;

% use fzero to find the rate where present value = 0
i_star = fzero(cf_sum, i_init);

% display result, %
fprintf('i_star: %.2f%%\n', i_star * 100);

% i_star is 19.17%
% this code is turning the cash flow into one large polynomial that is a function of r, standing in for ROR.
% the polynomial is then solved, with the result being the interest rate.