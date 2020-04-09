#####################################
## Audio Data Processing Functions ##
#####################################
library(tuneR)
library(signal)
## Convert frequency in Hz to MIDI pitch 0 - 127 (integer)
## 60 = Middle C in MIDI pitch
## A4 = reference frequency for A4, default is 440 Hz (69 is MIDI pitch of A4)
## Equal tempermant
## 100 cents between each MIDI pitch (integer), round to nearest whole number
freq2midi <- function(freq, A4 = 440){
  p <- round(69 + 12*log2(freq/A4))
  ## Check for 0
  if(any(p <= 0)){
    p[which(p <=0)] = 0
  }
  return(p)
}

## Convert MIDI pitch to frequency
midi2freq <- function(pitch, A4 = 440){
  f <- A4*2^((pitch - 69)/12)
  return(f)
}

## Convert MIDI pitch (0-127) to chroma number, 
## ref  is MIDI pitch of lowest note that chroma starts with, default is A0
# A0 = 21
# B0 = 23
# C1 = 24
# D1 = 26
# E1 = 28
# F1 = 29
# G1 = 31
midi2chroma <- function(midi, ref = 21){
  chrom <- (midi - ref)%%12 + 1
  labs <- c("A", "A#/Bb", "B", "C", "C#/Db", "D", "D#/Eb", "E", "F", 
            "F#/Gb", "G", "G#/Ab")
  ind <- ref - 21 + 1
  ## Reorder labels
  labels <- c(labs[ind:length(labs)], labs[1:ind-1])
  return(list(chroma = chrom, labels = labels))
}

## Test twinkle, twinkle, little star
#melody <- c(60, 60, 67, 67, 69, 69, 67, 65, 65, 64, 64, 62, 62, 60)
#all(freq2midi(midi2freq(melody)) == melody)
#midi2chroma(melody, ref = 24)

# Function to convert input audio to spectrogram, MIDI pitch representation and chromogram
## Audio = input audio signal
## sr = sampling rate of audio
## window.length = window duration in number of points for FFT, 
## size of the Fourier transform window
## step.time = step between windows in seconds for FFT
## A4 = frequency of pitch A4 in Hz, default is 440 Hz
## ref.pitch = lowest pitch on chromagram in MIDI pitch, default is 21 (A)
## plot = logical, display plots of spectrogram, MIDIgram, chromogram
encode <- function(audio, sr, window.length, A4 = 440, ref.pitch = 21, plot = FALSE){
  ## First find power spectrum
  #window.length <- sr*window.time/2
  spec <- specgram(audio, n = window.length, Fs = sr)
  # keep magnitude only of spectrogram (power)
  P <- abs(spec$S)
  
  
  ## Aggregate magnitude and frequencies by MIDI pitch value
  midi.gram <- rowsum(abs(spec$S), group = freq2midi(spec$f), reorder = TRUE)
  
  
  
  
  ## Aggregate magnitude by chroma values
  mc <- midi2chroma(as.numeric(row.names(midi.gram)), ref = ref.pitch)
  chroma <- mc$chroma
  labels <- mc$labels
  chromagram <- rowsum(midi.gram, group = chroma, reorder = TRUE)
  
  ## Normalize and convert to dB
  #P <- (P/max(P)
  P <- 10*log10(P+1)
  # normalize and convert to dB
  #midi.gram <- 10*log10(midi.gram + 1)
  #midi.gram <- midi.gram/max(midi.gram)
  #chromagram <- 10*log10(chromagram + 1)
  
  
  ## Plot
  if(plot){
    ## Spectrogram
  
    Spec.plot <- melt(t(P))
    Spec.plot$Freq <- rep(spec$f, each = length(spec$t))
    Spec.plot$Time <- rep(spec$t, times = length(spec$f))
    print(ggplot(Spec.plot, aes(Time, Freq, fill = value/max(value))) + geom_raster() +
      ggtitle('Spectrogram') + 
      scale_fill_continuous(type = "viridis") + 
      guides(fill=guide_legend(title="dB")) +
      ylab("Frequency [Hz]") + 
      xlab("Time [s]"))
    
    ## MIDI gram
    midi.plot <- melt(t(midi.gram))
    midi.plot$Freq <- rep(as.numeric(row.names(midi.gram)), each = length(spec$t))
    midi.plot$Time <- rep(spec$t, times = dim(midi.gram)[1])
    print(ggplot(midi.plot, aes(Time, Freq, fill = value/max(value))) + geom_raster() +
      ggtitle('MIDIgram') + 
      scale_fill_continuous(type = "viridis") + 
      guides(fill=guide_legend(title="dB")) + 
      scale_y_continuous(breaks = seq(0, 120, by = 10)) + 
      ylab("MIDI Pitch") + 
      xlab("Time [s]"))
    
    ## Chromagram
    labels <- factor(labels, levels = labels)
    chroma.plot <- melt(t(chromagram))
    chroma.plot$Freq <- rep(labels, each = length(spec$t))
    chroma.plot$Time <- rep(spec$t, times = dim(chromagram)[1])
    print(ggplot(chroma.plot, aes(Time, Freq, fill = value/max(value))) + geom_raster() +
      ggtitle('Chromagram') + 
      scale_fill_continuous(type = "viridis") + 
      guides(fill=guide_legend(title="dB")) + 
      #scale_y_discrete(breaks = seq(1, 12, by = 1)) + 
      ylab("Note Pitch") + 
      xlab("Time [s]"))
  }
  return(list(spec = spec, midi.gram = midi.gram, chromagram = chromagram, labels = labels))
}

## Plots

## plot a spectrogram, spec = spectrogram
plot.spec <- function(spec){
  P <- abs(spec$S)
  P <- 10*log10(P+1)
  Spec.plot <- melt(t(P))
  Spec.plot$Freq <- rep(spec$f, each = length(spec$t))
  Spec.plot$Time <- rep(spec$t, times = length(spec$f))
  print(ggplot(Spec.plot, aes(Time, Freq, fill = value/max(value))) + geom_raster() +
          ggtitle('Spectrogram') + 
          scale_fill_continuous(type = "viridis") + 
          guides(fill=guide_legend(title="dB")) +
          ylab("Frequency [Hz]") + 
          xlab("Time [s]"))
}


## MIDI gram
## tt = time points from spectrogram
plot.midi <- function(midi.gram, tt){
  midi.plot <- melt(t(midi.gram))
  midi.plot$Freq <- rep(as.numeric(row.names(midi.gram)), each = length(tt))
  midi.plot$Time <- rep(tt, times = dim(midi.gram)[1])
  print(ggplot(midi.plot, aes(Time, Freq, fill = value/max(value))) + geom_raster() +
          ggtitle('MIDIgram') + 
          scale_fill_continuous(type = "viridis") + 
          guides(fill=guide_legend(title="dB")) + 
          scale_y_continuous(breaks = seq(0, 120, by = 10)) + 
          ylab("MIDI Pitch") + 
          xlab("Time [s]"))
}


## Chromagram
## tt = time points from spectrogram
plot.chrom <- function(chromagram, labels, tt){
  labels <- factor(labels, levels = labels)
  chroma.plot <- melt(t(chromagram))
  chroma.plot$Freq <- rep(labels, each = length(tt))
  chroma.plot$Time <- rep(tt, times = dim(chromagram)[1])
  print(ggplot(chroma.plot, aes(Time, Freq, fill = value/max(value))) + geom_raster() +
          ggtitle('Chromagram') + 
          scale_fill_continuous(type = "viridis") + 
          guides(fill=guide_legend(title="dB")) + 
          #scale_y_discrete(breaks = seq(1, 12, by = 1)) + 
          ylab("Note Pitch") + 
          xlab("Time [s]"))
}
 

## Calculate spectral flatness by frame, can't be on log10 scale, 
#uses sfm from seewave package
spec.flatness <- function(spectrogram){
  ## Spectrogram should be magnitude, not on log scale (no dB)
  ## Spectrogram is Hz x # frames
  sf <- list()
  for(i in 1:dim(spectrogram)[2]){
    sf[i] <- sfm(spectrogram[,i])
  }
  return(unlist(sf))
}



###########################################
## Process Beethoven and Produce Curves ##
##########################################


## Main function to calculate tempo, volume and spectral flatness curves
## Inputs: file.list is a list of audio filenames to read in and process
## Assumes each filename corresponds to ONE complete movement for ONE orchestra
## folder is the folder path to the files in file.list
## In order of file.list order (i.e. alphabetical)
## A4.list is a list of the tuning pitches of A4 for each orchestra in orchs
## ref.word is which string in each filename should be treated as the "reference" recording
## to align the other recordings to
## Assumes reference recording contains "MIDI" in the filename
## f.inds1 is list of length = # of files of indices of beginning of each frequency group
## f.inds2 is list of length = # of files of indices of ending of each frequency group
## can change frequencies of f.inds1 and f.inds2 to account for different tunings
## of orchestras
## nMFCC is number of MFCC coefficients to calculate
## wintime is window length in s for MFCC
## hoptime is steps in s between windows for MFCC
## window.length is the size of the Fourier transform window
## 4410*2 gives a frequency resolution of 5 Hz and time resolution of 0.1 s
## Returns: deriv.df a list of the tempo, volume, spectral flatness, volume by frequency,
## spectral flatness by frequency and MFCC curves
## and the tempo alignment indices for the reference recording (align.ref)
## and alignment indices for each orchestra (align.orch)
## Run once for each movement/piece

process.curves <- function(file.list, folder, A4.list, ref.word = "MIDI",
                           f.inds1, f.inds2, nMFCC, wintime = 0.05, hoptime = 0.1,
                           window.length = 4410*2){
  chrom <- list() ## chromagram
  sf <- list() ## spectral flatness
  mfcc.list <- list() ## mfcc
  spec.list <- list() ## store spectrograms
  count <- 1
  for(i in file.list){
    print(i)
    filename <- paste(folder, i, sep = '') ## read in each audio file
    curr <- readMP3(filename)
    sr <- curr@samp.rate
    pieces <- curr@left
    ind <- which(pieces > 0 )
    ind1 <- ind[1] ## remove silence before piece
    ind2 <- tail(ind, n = 1) ## remove silence from end of piece
    curr <- pieces[ind1:ind2]
    e <- encode(curr, sr, window.length, A4 = A4.list[count])
    spec.list[[count]] <- abs(e$spec$S)
    chrom[[count]] <- e$chromagram # store chromagram
    sf[[count]] <- spec.flatness(abs(e$spec$S)) # calculate and store spectral flatness
    ## MFCC - all files
    mc <- melfcc(readMP3(filename), sr = sr, wintime = wintime, hoptime = hoptime)
    mc[which(is.nan(mc))] <- 0 ##set Nan to 0
    mfcc.list[[count]] <- mc
    count <- count + 1
  
  }
  
  ## Find which filename should be treated as reference recording
  ref.ind <- which(grepl(ref.word, file.list))
  query.ind <- which(c(1,2) != ref.ind)
  
  
  # (1) Align and store indices
  align <- dtw(t(chrom[[ref.ind]]), t(chrom[[query.ind]]), 
               step.pattern = asymmetric, dist.method = "cosine")
  if(length(align$index2) <= 500){
    n.smooth <- 71
  } else{
    n.smooth <- 501
  }
  
  ind.df <- align$index2
  ref.ind.df <- align$index1
  
  
  # (2) Store Tempo Curves (find derivative with Savitzky-Golay filter)
  dt <- sgolayfilt(align$index2, m = 1, n = n.smooth)
  
  ## Normalize area under curve
  x <- align$index1/max(align$index1)
  y <- dt
  AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
  
  tempo.df <- dt/AUC
  
  # (3) Volume Curve - Overall
  vol.chrom <- chrom[[query.ind]][, align$index2] ## align
  avg <- mean(colMeans(vol.chrom))
  vol.curve <- colMeans(vol.chrom)/avg # Normalize by average volume
  
  x <- align$index1/max(align$index1)
  y <- vol.curve
  AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
  
  volume.df <- vol.curve/AUC
  
  
  # (4) Spectral Flatness - Just Subset and Normalize
  sf.chrom <- sf[[query.ind]][align$index2] ## align

  x <- align$index1/max(align$index1)
  y <- sf.chrom
  AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
  
  sf.df <- sf.chrom/AUC
  
  # (5) Volume by Frequency Group 
  # (6) Spectral Flatness by Frequency Group
  ## Matrix to save frequency curves
  vol.chrom <- spec.list[[query.ind]][, align$index2] ## align
  avg <- mean(colMeans(vol.chrom))
  avg.sf <- mean(spec.flatness(vol.chrom))
  x <- align$index1/max(align$index1)
  vol.freq <- matrix(0, nrow = length(f.inds1), ncol = ncol(vol.chrom))
  sf.freq <- matrix(0, nrow = length(f.inds1), ncol = ncol(vol.chrom))
  for(j in 1:length(f.inds1)){
    vol.curve <- colMeans(vol.chrom[f.inds1[j]:f.inds2[j], ])/avg
    y <- vol.curve
    AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
    vol.freq[j, ] <- vol.curve/AUC
    
    
    sf.curve <- spec.flatness(vol.chrom[f.inds1[j]:f.inds2[j], ])/avg.sf
    y <- sf.curve
    AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
    sf.freq[j, ] <- sf.curve/AUC
  }
  
  # (7) MFCC
  align2 <- dtw(mfcc.list[[ref.ind]], mfcc.list[[query.ind]], 
               step.pattern = asymmetric, dist.method = "cosine")
  m.chrom <- t(mfcc.list[[query.ind]])[, align2$index2] ## align
  for(m in 1:(nMFCC-1)){
    x <- align2$index1/max(align2$index1)
    y <- m.chrom[m, ]
    AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
    m.chrom[m,] <- y/AUC
  }
  ## Store MFCCs
  mfcc.df <- m.chrom
  
  
  rm(spec.list, chrom, sf, mfcc.list)
  
  
  deriv.df <- list(tempo = tempo.df, volume = volume.df,
                   SF = sf.df, vol.freq = vol.freq, sf.freq = sf.freq, mfcc = m.chrom,
                   align.ref = ref.ind.df, align.orch = ind.df)
  
  return(deriv.df)
}





## Plot spectrogram curves

## chrom.list is a nested list, first level is each orchestra, named
## Second level of list is the movement for each symphony
## Values in the list are matrices of volume curves by frequency
## r is row of volume matrix/frequency group to select
## Returns list of length = number of movements of plots of the curve in chrom.list
## i.e. if chrom.list contains tempo curve values, then plot.list is a length 4
## list of plots of the tempo curves for each orchestra by each movement
plot.spec.curves <- function(chrom.list, r){
  plot.list <- list()
  for(i in 1:length(chrom.list[[1]])){
    m1 <- data.frame(lapply(lapply(chrom.list, `[[`, i), `[`, r,))
    n <- nrow(m1)
    b.m1 <- melt(m1)
    colnames(b.m1) <- c("Orchestra", "Density")
    b.m1$Index <- rep(seq(n), times = length(orchs))/n
    p <- ggplot(b.m1, aes(x = Index, y = Density, col = Orchestra)) +
      geom_line() +
      ggtitle(paste("Mvmt. ", i))
    plot.list[[i]] <- p
  }
  return(plot.list)
  
}



#############
### Timbre ##
#############


## Process spectrograms, too large to store multiple in memory at the same time
## So align chromagrams, then using aligned indices align and subset spectrograms
## Group spectrograms by frequency and sum to form volume curves by frequency group
## Captures aspects of instrument balance and timbre
## Inputs: file.list is a list of audio filenames to read in and process
## reference recording not included in file.list
## Assumes each filename corresponds to ONE complete movement for ONE orchestra
## folder is the folder path to the files in file.list
## orchs is a list of the different orchestras contained in file.list to compare
## In order of file.list order (i.e. alphabetical)
## A4.list is a list of the tuning pitches of A4 for each orchestra in orchs
## length(orchs) == length(A4.list)
## ref.word is which string in each filename should be treated as the "reference" recording
## window.length is the size of the Fourier transform window
## 4410*2 gives a frequency resolution of 5 Hz and time resolution of 0.1 s
## f.inds1 is list of length = # of files of indices of beginning of each frequency group
## f.inds2 is list of length = # of files of indices of ending of each frequency group
## can change frequencies of f.inds1 and f.inds2 to account for different tunings
## of orchestras
## align is a list (of same length as file.list and in same order) 
## of alignment indices for each orchestra/movement
## nMove is the number of movements
## nMFCC is number of MFCC coefficients to calculate
## wintime is window length in s for MFCC
## hoptime is steps in s between windows for MFCC
## Returns: sf.df a nested list of matrices
## First layer of list is each orchestra, then one matrix for each orchestra
## Dimensions are number of frequency groups (length f.inds1) x number of time points 
## mfcc.df is same nested list of frequencies, nMFCC-1 x number of time points
## First MFCC is dropped
process.timbre <- function(file.list, folder, orchs, A4.list, ref.word,
                         window.length = 4410*2, f.inds1, f.inds2, align, nMove,
                         nMFCC, wintime = 0.05, hoptime = 0.1){
  chrom <- list() ## Store spectral flatness by frequency group
  mfcc.list <- list() ## Store MFCC coefficients
  count <- 1
  count.mfcc <- 1
  for(i in file.list){
    print(i)
    filename <- paste(folder, i, sep = '')
    
    if(!grepl(ref.word, filename)){ ## Only need SF for non-MIDI files
      curr <- readMP3(filename)
      sr <- curr@samp.rate
      pieces <- curr@left
      ind <- which(pieces > 0 )
      ind1 <- ind[1] ## remove silence before piece
      ind2 <- tail(ind, n = 1) ## remove silence from end of piece
      curr <- pieces[ind1:ind2]
      e <- encode(curr, sr, window.length, A4 = A4.list[count])
      
      ## Matrix to save frequency curves
      vol.chrom <- abs(e$spec$S)[, align[[count]]] ## align temporally
      avg <- mean(spec.flatness(vol.chrom))
      x <- align[[count]]/max(align[[count]])
      chrom[[count]] <- matrix(0, nrow = length(f.inds1[[1]]), ncol = ncol(vol.chrom))
      for(j in 1:length(f.inds1[[count]])){
        sf.curve <- spec.flatness(vol.chrom[f.inds1[[count]][j]:f.inds2[[count]][j], ])/avg
        
        y <- sf.curve
        AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
        chrom[[count]][j, ] <- sf.curve/AUC
      }
      count <- count + 1
    } 
    
    ## MFCC - all files
    mc <- melfcc(readMP3(filename), sr = sr, wintime = wintime, hoptime = hoptime)
    mc[which(is.nan(mc))] <- 0 ##set Nan to 0
    mfcc.list[[count.mfcc]] <- mc
    count.mfcc <- count.mfcc + 1

  }
  
  
  sf.df <- vector("list", length(orchs))
  count <- 1
  new.list <- file.list[-grep(ref.word, file.list)]
  for(orc in orchs){
    curr.orc <- grep(orc, new.list) ## Drop MIDI from file.list to index correctly
    for(i in 1:nMove){
      j <- curr.orc[i]
      sf.df[[count]][[i]] <- chrom[[j]]
    }
    count <- count + 1
  }
  
  ## Find which filenames should be treated as reference recordings
  ref.ind <- which(grepl(ref.word, file.list))
  ## Align and store curves of interest
  mfcc.df <- vector("list", length(orchs))
  count <- 1
  for(orc in orchs){
    curr.orc <- grep(orc, file.list)
    for(i in 1:length(ref.ind)){
      
      j <- curr.orc[i]
      # (1) Align 
      align <- dtw(mfcc.list[[ref.ind[i]]], mfcc.list[[j]], 
                   step.pattern = asymmetric, dist.method = "cosine")
      # (4) Align and Normalize MFCC
      m.chrom <- t(mfcc.list[[j]])[, align$index2] ## align
      for(m in 1:(nMFCC-1)){
        x <- align$index1/max(align$index1)
        y <- m.chrom[m, ]
        AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
        m.chrom[m,] <- y/AUC
      }
      ## Store MFCCs
      mfcc.df[[count]][[i]] <- m.chrom
    }
    count <- count + 1
    print(curr.orc)
  }
  
  rm(chrom, mfcc.list)
  sf.df <- setNames(sf.df, orchs)
  mfcc.df <- setNames(mfcc.df, orchs)
  return(list(SF.byF = sf.df, MFCC = mfcc.df))
}

## Process spectrograms, too large to store multiple in memory at the same time
## So align chromagrams, then using aligned indices align and subset spectrograms
## Group spectrograms by frequency and sum to form volume curves by frequency group
## Captures aspects of instrument balance and timbre
## Inputs: file.list is a list of audio filenames to read in and process
## reference recording not included in file.list
## Assumes each filename corresponds to ONE complete movement for ONE orchestra
## folder is the folder path to the files in file.list
## orchs is a list of the different orchestras contained in file.list to compare
## In order of file.list order (i.e. alphabetical)
## A4.list is a list of the tuning pitches of A4 for each orchestra in orchs
## length(orchs) == length(A4.list)
## ref.word is which string in each filename should be treated as the "reference" recording
## window.length is the size of the Fourier transform window
## 4410*2 gives a frequency resolution of 5 Hz and time resolution of 0.1 s
## f.inds1 is list of length = # of files of indices of beginning of each frequency group
## f.inds2 is list of length = # of files of indices of ending of each frequency group
## can change frequencies of f.inds1 and f.inds2 to account for different tunings
## of orchestras
## align is a list (of same length as file.list and in same order) 
## of alignment indices for each orchestra/movement
## nMove is the number of movements
## nMFCC is number of MFCC coefficients to calculate
## wintime is window length in s for MFCC
## hoptime is steps in s between windows for MFCC
## Returns: sf.df a nested list of matrices
## First layer of list is each orchestra, then one matrix for each orchestra
## Dimensions are number of frequency groups (length f.inds1) x number of time points 
## mfcc.df is same nested list of frequencies, nMFCC-1 x number of time points
## First MFCC is dropped
process.timbre9 <- function(file.list, folder, orchs, A4.list, ref.word,
                           window.length = 4410*2, f.inds1, f.inds2, align, nMove,
                           nMFCC, wintime = 0.05, hoptime = 0.1){
  chrom <- list() ## Store spectral flatness by frequency group
  mfcc.list <- list() ## Store MFCC coefficients

  count <- 1
  count.mfcc <- 1
  for(i in file.list){
    if(length(i) == 1){ ## recording in file.list is complete movement
      print(i)
      filename <- paste(folder, i, sep = '')
      if(!grepl(ref.word, filename)){ ## Only need SF for non-MIDI files
        curr <- readMP3(filename)
        sr <- curr@samp.rate
        pieces <- curr@left
        ind <- which(pieces > 0 )
        ind1 <- ind[1] ## remove silence before piece
        ind2 <- tail(ind, n = 1) ## remove silence from end of piece
        curr <- pieces[ind1:ind2]
        e <- encode(curr, sr, window.length, A4 = A4.list[count])
        vol <- abs(e$spec$S)
        
        ## Matrix to save frequency curves
        ## Matrix to save frequency curves
        vol.chrom <- vol[, align[[count]]] ## align temporally
        avg <- mean(spec.flatness(vol.chrom))
        x <- align[[count]]/max(align[[count]])
        chrom[[count]] <- matrix(0, nrow = length(f.inds1[[1]]), ncol = ncol(vol.chrom))
        for(j in 1:length(f.inds1[[count]])){
          sf.curve <- spec.flatness(vol.chrom[f.inds1[[count]][j]:f.inds2[[count]][j], ])/avg
          
          y <- sf.curve
          AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
          chrom[[count]][j, ] <- sf.curve/AUC
        }
        count <- count + 1
      }
      ## MFCC - all files
      mc <- melfcc(readMP3(filename), sr = sr, wintime = wintime, hoptime = hoptime)
      mc[which(is.nan(mc))] <- 0 ##set Nan to 0
      mfcc.list[[count.mfcc]] <- mc
      
      
    } else{ ## need to combine movements
      move4 <- list()
      move4.mfcc <- list()
      cc <- 1
      for(j in i){
        print(j)
        filename <- paste(folder, j, sep = '')
        #curr <- readMidi(filename)
        curr <- readMP3(filename)
        sr <- curr@samp.rate
        pieces <- curr@left
        ind <- which(pieces > 0 )
        ind1 <- ind[1] ## remove silence before piece
        ind2 <- tail(ind, n = 1) ## remove silence from end of piece
        curr <- pieces[ind1:ind2]
        e <- encode(curr, sr, window.length, A4 = A4.list[count])
        move4[[cc]] <- abs(e$spec$S)
        
        ## MFCC
        mc <- melfcc(readMP3(filename), sr = sr, wintime = wintime, hoptime = hoptime)
        mc[which(is.nan(mc))] <- 0 ##set Nan to 0
        move4.mfcc[[cc]] <- mc
        
        cc <- cc + 1
      }
      vol <- do.call(cbind, move4) ## Concatenate
      mfcc.list[[count.mfcc]] <- do.call(rbind, move4.mfcc) ## Concatenate

      ## Matrix to save frequency curves
      vol.chrom <- vol[, align[[count]]] ## align temporally
      avg <- mean(spec.flatness(vol.chrom))
      x <- align[[count]]/max(align[[count]])
      chrom[[count]] <- matrix(0, nrow = length(f.inds1[[1]]), ncol = ncol(vol.chrom))
      for(j in 1:length(f.inds1[[count]])){
        sf.curve <- spec.flatness(vol.chrom[f.inds1[[count]][j]:f.inds2[[count]][j], ])/avg
        
        y <- sf.curve
        AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
        chrom[[count]][j, ] <- sf.curve/AUC
      }
      count <- count + 1
    }
    
    count.mfcc <- count.mfcc + 1
  }
  
  
  
  
  sf.df <- vector("list", length(orchs))
  count <- 1
  new.list <- file.list[-grep(ref.word, file.list)]
  for(orc in orchs){
    curr.orc <- grep(orc, new.list) ## Drop MIDI from file.list to index correctly
    for(i in 1:nMove){
      j <- curr.orc[i]
      sf.df[[count]][[i]] <- chrom[[j]]
    }
    count <- count + 1
  }
  
  ## Find which filenames should be treated as reference recordings
  ref.ind <- which(grepl(ref.word, file.list))
  ## Align and store curves of interest
  mfcc.df <- vector("list", length(orchs))
  count <- 1
  for(orc in orchs){
    curr.orc <- grep(orc, file.list)
    for(i in 1:length(ref.ind)){
      
      j <- curr.orc[i]
      # (1) Align 
      align <- dtw(mfcc.list[[ref.ind[i]]], mfcc.list[[j]], 
                   step.pattern = asymmetric, dist.method = "cosine")
      # (4) Align and Normalize MFCC
      m.chrom <- t(mfcc.list[[j]])[, align$index2] ## align
      for(m in 1:(nMFCC-1)){
        x <- align$index1/max(align$index1)
        y <- m.chrom[m, ]
        AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
        m.chrom[m,] <- y/AUC
      }
      ## Store MFCCs
      mfcc.df[[count]][[i]] <- m.chrom
    }
    count <- count + 1
    print(curr.orc)
  }
  
  rm(chrom, mfcc.list)
  sf.df <- setNames(sf.df, orchs)
  mfcc.df <- setNames(mfcc.df, orchs)
  return(list(SF.byF = sf.df, MFCC = mfcc.df))
}
