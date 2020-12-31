library(bitops)

##The algorithm is followed from: Calculation of Mappings Between
##One and n-dimensional Values Using the Hilbert Space-filling Curve, 
##J.K. Lawder. The functions correspond to the algorithm step in paper.

##Key function: convert
##Key parameters: input: the point (1-d) needs to tranform; 
##order: interger, the order of Hilbert Curve, must be less than 31 cause 
##       the bit limit of R; 
##dim: integer, the dim after transformation (the dimension of Hilbert Curve);
##random: logic, whether a random setting need, for RQMC setting only.

DerivedKey = function(input, order, dim)
{
  ORDER = max(order, 1)
  tmp = bitShiftL(1, dim)
  rho = rep(0, order)
  
  for (i in 1:ORDER)
  {
    input = input*tmp
    rho[i] = as.integer(input)
    input = input - rho[i]
  }
  
  return(rho)
}

calc_J = function(P, dim)
{
  J = dim
  k = 0
  for (i in 1:(dim-1))
  {
    k = i
    if ((bitAnd(bitShiftR(P, i), 1)) == (bitAnd(P, 1)))
    {
      next
    }
    else
    {
      break
    }
  }
  
  if (i != dim)
  {
    J = J - i
  }
  return(J)
}

calc_T = function(P)
{
  if (P < 3)
  {
    return(0)
  }
  
  if (P%%2 == 0)
  {
    return(bitXor(P - 1, (P - 1)/2))
  }
  
  return(bitXor(P - 2, (P - 2)/2))
}

calc_tS_tT = function(xJ, val, dim)
{
  retval = val
  if (xJ%%dim != 0)
  {
    temp1 = bitShiftR(val, xJ%%dim)
    temp2 = bitShiftL(val, dim-(xJ%%dim))
    retval = bitOr(temp1, temp2)
    retval = bitAnd(retval, bitShiftL(1, dim) - 1)
  }
  
  return(retval)
}

convert = function(input, order, dim, random = F)
{
  ORDER = max(1, order)
  rho = DerivedKey(input, order, dim)
  
  xJ = 0
  W = 0
  pt = rep(0, dim)
  
  mask = 0.5
  for (i in 1:ORDER)
  {
    P = rho[i]
    S = bitXor(P, P/2)
    tS = calc_tS_tT(xJ, S, dim)
    A = bitXor(W, tS)
    
    j = dim
    while(j >= 1 && A > 0)
    {
      if (bitAnd(A, 1))
      {
        pt[j] = pt[j] + mask
      }
      j = j - 1
      A = bitShiftR(A, 1)
    }
    
    T1 = calc_T(P)
    tT = calc_tS_tT(xJ, T1, dim)
    W = bitXor(W, tT)
    J = calc_J(P, dim)
    xJ = xJ + J - 1
    mask = mask/2
  }
  
  if (random == T)
  {
    for (b in 1:dim)
    {
      pt[b] = pt[b] + runif(1)/(bitShiftL(1, order))
    }
  }
  
  return(pt)
}


