# Complete-Root-Locus-Tool-v1 
# An enhanced Python root locus plotting program featuring:
  - open loop poles and zeros
  - complete root locus plot for both k > 0 and k < 0
  - angles and intersection of asymptotes
  - candidates of breakaway points solved by d/ds(loop gain) = 0
  - marginally stable poles location with corresponding k's
# Demonstration
![Complete_Root_Locus_tool_v1_2024_1126_011348](https://github.com/user-attachments/assets/4d15991c-911e-4911-bfb0-8e96e8da6499)

# User guide
  - Transfer Function format:
    - forward transfer function       = G(s)  = GN(s)/GD(s)
    - feedback transfer function      = H(s)  = HN(s)/HD(s)
    - loop gain                       = GL(s) = G(s) * H(s) 

  - Entering Coefficients
    - Please enter the 1 or 2 degree factors in the format of coefficients in the corresponding csv files.
    - eg. s(s+2)(s^2+3)
      - csv format:
        -    0     1   0   (s)
        -    0     1   2   (s+2)
        -    1     0   3   (s^2 + 3)

  - csv filename matching list:
    - GN(s) => G(s)_nominator.csv
    - GD(s) => G(s)_denominator.csv
    - HN(s) => H(s)_nominator.csv
    - HD(s) => H(s)_denominator.csv

  - Preference settings: settings.csv
    - parameters:
    - prec     = displayed decimal precision
    - xRange   = displayed x-axis value range in both +/- direcitons. Set to -1 for automatically adjusting
      - y axis will fit automatically with same unit length under all circumstances
    - lb       = k(gain) lowerbound = k0 * 10^lb		
    - ub       = k(gain) upperbound = k0 * 10^ub
      - k0 = the gain such that the highest deg terms have equal coefficient in both GL(s)'s nominator and denominator
    - quantity = sample quantity
