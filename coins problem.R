# coins problem

PDcoins_raw <- function(Coins, MaxValue) {
  ValueArray <- rep(Inf, MaxValue + 1)  # Initialize with large numbers
  ValueArray[1] <- 0  # Base case: 0 coins needed for value 0
  for (Coin in Coins) {
    if (Coin <= MaxValue) {
      ValueArray[Coin] <- 1
    }
  }
  for (Value in 1:MaxValue) {
    if (ValueArray[Value] < Inf) {  # Only proceed if the value is reachable
      for (SumCoin in Coins) {
        NewValue <- Value + SumCoin
        if (NewValue <= MaxValue) {
          ValueArray[NewValue] <- min(ValueArray[NewValue], ValueArray[Value] + 1)
        }
      }
    }
  }
  
  return(ifelse(ValueArray[MaxValue] == Inf, -1, ValueArray[MaxValue])) # -1 means unreachable
}

Coins <- c(3,2)
MaxValue <- 24
NumCoins <- PDcoins_raw(Coins, MaxValue)
print(NumCoins)


PDcoins <- function(Coins, MaxValue) {
  CoinCombination <- vector("list", MaxValue)  # Initialize as list
  ValueArray <- rep(Inf, MaxValue + 1)  # Initialize with large numbers  # Base case: 0 coins needed for value 0
  for (Coin in Coins) {
    if (Coin <= MaxValue) {
      ValueArray[Coin] <- 1
    }
  }
  # Pre-populate worst-case scenario for each value: maximum number of 1's
  MinCoin <- min(Coins)
  for (m in 1:MaxValue) {
    CoinCombination[[m]] <- rep(MinCoin, m)  # Assume worst case is using all 1's
  }
  
  for (Coin in Coins) {
    if (Coin <= MaxValue) {
      ValueArray[Coin] <- 1  # Directly assign single-coin cases
      CoinCombination[[Coin]] <- c(Coin)  # If a value matches a coin exactly, use only that coin
    }
  }
  for (Value in 1:MaxValue) {
    if (ValueArray[Value] < Inf) {  # Only proceed if the value is reachable
      for (SumCoin in Coins) {
        NewValue <- Value + SumCoin
        if (NewValue <= MaxValue) {
          oldCombination <- as.numeric(c(unlist(CoinCombination[[Value]])))
          NewCombination <- c(oldCombination , SumCoin)
          
          # Update only if the new combination is shorter than the current
          if (length(NewCombination) < length(CoinCombination[[NewValue]])) {
            ValueArray[NewValue] <- ValueArray[Value] + 1
            CoinCombination[[NewValue]] <- NewCombination  # Keep the shortest valid combination
          }
        }
      }
    }
  }
  ValueArray[MaxValue] <- ifelse(ValueArray[MaxValue] == Inf, -1, ValueArray[MaxValue])
  if (length(CoinCombination[[MaxValue]]) == MaxValue) {
    CoinCombination[[MaxValue]] <- -1
  }
  return(list(MinNumCoins = ValueArray[MaxValue], CoinCombination = CoinCombination[[MaxValue]]))
  
}

Coins <- c(2,3)
MaxValue <- 24
time_taken <- system.time(Result <- PDcoins(Coins, MaxValue))
print(Result)
print(time_taken)
