# Complete-Root-Locus-Tool-v1 
# An enhanced Python root locus plotting program featuring:
  - open loop poles and zeros
  - complete root locus plot for both k > 0 and k < 0
  - angles and intersection of asymptotes
  - candidates of breakaway points solved by d/ds(loop gain) = 0
  - marginally stable poles location with corresponding k's
# Plot Demonstration
![Complete_Root_Locus_tool_v1_2024_1126_011348](https://github.com/user-attachments/assets/4d15991c-911e-4911-bfb0-8e96e8da6499)

# User guide
  - Software
    - Complete Root Locus tool.py (main)
    - CompleteRL.py (calculation module)

  - csv files :
    - Transfer function G(s) and H(s)
      - G(s)_nominator.csv
        - GN(s)
      - G(s)_denominator.csv
        - GD(s)
      - H(s)_nominator.csv
        - HN(s)
      - H(s)_denominator.csv
        - HD(s)
    - Preference Settings
      - settings.csv
        
  - Transfer Function format:
    - forward transfer function       = G(s)  = GN(s)/GD(s)
    - feedback transfer function      = H(s)  = HN(s)/HD(s)
    - loop gain                       = GL(s) = G(s) * H(s) 

    - Entering Coefficients
      - Please enter the 1 or 2 degree factors in the format of coefficients in the corresponding csv files.
      - eg. s(s+2)(s^2+3)
        - entering format in csv file:
          -    0     1   0   (s)
          -    0     1   2   (s+2)
          -    1     0   3   (s^2 + 3)

  - Preference settings: settings.csv
    - parameters:
    - prec     = displayed decimal precision
    - xRange   = displayed x-axis value range in both +/- direcitons. Set to -1 for automatically adjusting
      - y axis will fit automatically with same unit length under all circumstances
    - lb       = k(gain) lowerbound = k0 * 10^lb		
    - ub       = k(gain) upperbound = k0 * 10^ub
      - k0 = the gain such that the highest deg terms have equal coefficient in both GL(s)'s nominator and denominator
    - quantity = sample quantity
