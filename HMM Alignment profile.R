### HMM Alignment profile 

read_profile_input <- function(path) {
  lines <- trimws(readLines(path))
  # first non-empty line is threshold
  lines <- lines[nzchar(lines)]
  x <- lines[[1]]
  theta <- as.numeric(unlist(strsplit(lines[sep_idx[1] + 1]," ")))
  pseudocounts <- theta[[2]]
  theta <- theta[[1]]
  # find separators '--------'
  sep_idx <- which(lines == "--------")
  if(length(sep_idx) < 2) stop("Unexpected input format: need two '--------' lines")
  alphabet <- strsplit(lines[sep_idx[2] + 1], "\\s+")[[1]]
  seq_profile <- lines[(sep_idx[3] + 1):length(lines)]
  seq_profile <- seq_profile[nzchar(seq_profile)]
  list(threshold = theta, alphabet = alphabet, seq_profile = seq_profile, pseudocounts = pseudocounts,x = x)
}

viterbi_profile <- function(hmm, x) {
  # hmm: list(transition = Tmat, emission = Emat)
  # x: character string (e.g., "AEFDFDC")
  states <- rownames(hmm$transition)
  trans <- hmm$transition
  emit  <- hmm$emission
  Lx <- nchar(x)
  seq_chars <- strsplit(x, "")[[1]]
  
  # identify emitting vs non-emitting states
  is_emitting <- rownames(emit) %in% states & rowSums(emit) > 0
  # ensure names align
  names(is_emitting) <- states
  
  # log matrices (use -Inf for zero)
  log_trans <- matrix(-Inf, nrow = nrow(trans), ncol = ncol(trans),
                      dimnames = dimnames(trans))
  nz <- trans > 0
  log_trans[nz] <- log(trans[nz])
  log_emit <- matrix(-Inf, nrow = nrow(emit), ncol = ncol(emit),
                     dimnames = dimnames(emit))
  nzE <- emit > 0
  log_emit[nzE] <- log(emit[nzE])
  
  # DP arrays: log-probabilities and backpointers
  # positions 0..Lx (0 = before consuming any char)
  nstates <- length(states)
  V <- matrix(-Inf, nrow = nstates, ncol = Lx + 1, dimnames = list(states, 0:Lx))
  prev_state <- matrix(NA_character_, nrow = nstates, ncol = Lx + 1, dimnames = list(states, 0:Lx))
  prev_pos   <- matrix(NA_integer_, nrow = nstates, ncol = Lx + 1, dimnames = list(states, 0:Lx))
  
  # helper to index state name -> row index
  state_idx <- function(s) which(states == s)
  
  # initialize: start at S at pos 0 with prob 1 (log 0)
  if(!("S" %in% states)) stop("HMM has no S state")
  V["S", "0"] <- 0
  
  # closure over silent transitions at a given pos:
  silent_closure <- function(pos) {
    changed <- TRUE
    while(changed) {
      changed <- FALSE
      for(u in states) {
        iu <- state_idx(u)
        if(is.infinite(V[iu, as.character(pos)])) next
        for(v in states) {
          # only consider silent (non-emitting) target states
          if(is_emitting[v]) next
          jv <- state_idx(v)
          lt <- log_trans[iu, jv]
          if(is.infinite(lt)) next
          cand <- V[iu, as.character(pos)] + lt
          if(cand > V[jv, as.character(pos)]) {
            V[jv, as.character(pos)] <<- cand
            prev_state[jv, as.character(pos)] <<- u
            prev_pos[jv, as.character(pos)] <<- pos
            changed <- TRUE
          }
        }
      }
    }
  }
  
  # initial silent closure at pos 0
  silent_closure(0)
  
  # main loop: consume characters
  for(pos in 0:(Lx-1)) {
    ch <- seq_chars[pos + 1]
    # for each emitting state e, compute best coming from any state u at pos
    for(e in states[is_emitting]) {
      je <- state_idx(e)
      best_val <- -Inf
      best_prev <- NA_character_
      best_prev_pos <- NA_integer_
      # emission probability for this char
      le <- log_emit[je, ch]
      if(is.infinite(le)) {
        # impossible emission -> skip (remains -Inf)
        next
      }
      for(u in states) {
        iu <- state_idx(u)
        if(is.infinite(V[iu, as.character(pos)])) next
        lt <- log_trans[iu, je]
        if(is.infinite(lt)) next
        cand <- V[iu, as.character(pos)] + lt + le
        if(cand > best_val) {
          best_val <- cand
          best_prev <- u
          best_prev_pos <- pos
        }
      }
      if(!is.infinite(best_val)) {
        V[je, as.character(pos + 1)] <- best_val
        prev_state[je, as.character(pos + 1)] <- best_prev
        prev_pos[je, as.character(pos + 1)] <- best_prev_pos
      }
    }
    # after filling emitting states at pos+1, run silent closure at pos+1
    silent_closure(pos + 1)
  }
  
  # finally, transition to End 'E' (may be silent). choose best u at pos Lx with trans[u,E]
  if(!("E" %in% states)) stop("HMM has no E state")
  best_val <- -Inf
  best_prev <- NA_character_
  best_prev_pos <- NA_integer_
  for(u in states) {
    iu <- state_idx(u)
    if(is.infinite(V[iu, as.character(Lx)])) next
    lt <- log_trans[iu, "E"]
    if(is.infinite(lt)) next
    cand <- V[iu, as.character(Lx)] + lt
    if(cand > best_val) {
      best_val <- cand
      best_prev <- u
      best_prev_pos <- Lx
    }
  }
  # if no path to E found, try to allow E itself if it had non -Inf at pos Lx (rare)
  if(is.infinite(best_val) && !is.infinite(V["E", as.character(Lx)])) {
    best_val <- V["E", as.character(Lx)]
    best_prev <- prev_state["E", as.character(Lx)]
    best_prev_pos <- prev_pos["E", as.character(Lx)]
  }
  if(is.infinite(best_val)) {
    stop("No valid path to End found for the given sequence")
  }
  
  # set E at pos Lx
  V["E", as.character(Lx)] <- best_val
  prev_state["E", as.character(Lx)] <- best_prev
  prev_pos["E", as.character(Lx)] <- best_prev_pos
  
  # backtrack from E at pos Lx to S at pos 0
  path_states <- character(0)
  cur_state <- "E"
  cur_pos <- Lx
  while(TRUE) {
    ps <- prev_state[state_idx(cur_state), as.character(cur_pos)]
    pp <- prev_pos[state_idx(cur_state), as.character(cur_pos)]
    # append current state (we will reverse later). Do not append S sentinel at the end.
    path_states <- c(path_states, cur_state)
    if(is.na(ps) || is.na(pp)) break
    cur_state <- ps
    cur_pos <- pp
    if(cur_state == "S") {
      path_states <- c(path_states, "S")
      break
    }
  }
  # reverse and drop S and E
  path_states <- rev(path_states)
  # remove leading S and trailing E if present
  if(length(path_states) > 0 && path_states[1] == "S") path_states <- path_states[-1]
  if(length(path_states) > 0 && tail(path_states, 1) == "E") path_states <- head(path_states, -1)
  
  # collapse consecutive identical states? No — keep full path including deletes and inserts.
  return(path_states)
}


### test

input_file <- choose.files()
input_lines <- read_profile_input(input_file)
hmm <- build_profile_hmm(input_lines$threshold,input_lines$alphabet,input_lines$seq_profile,
                         input_lines$pseudocounts)
optimal_path <- viterbi_profile(hmm,input_lines$x)
cat(optimal_path," ")
