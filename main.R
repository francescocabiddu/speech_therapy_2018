# load libraries
lib <- c("magrittr", "tidyverse", 
         "data.table", "fastmatch",
         "beepr", "mailR", "janitor")
lapply(lib, require, character.only = TRUE)
rm(lib)

# dfs are not updated in wordspace

#### convert reference dfs ####
# convert mot_tokens and chi_tokens "UU" to "AH" and recalculate phoneme frames
mot_tokens_converted <- mot_tokens %>%
  mutate(phon = phon %>%
           str_replace_all("UU", "AH")) %>%
  select(phon) %>%
  (function(y) {
    reference <- y[!duplicated(y$phon),] %>%
      mutate(split = str_split(phon, "_") %>%
               sapply(function(x) {
                 adj_comb_complete(x)
               })) %>%
      mutate(seq2 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 2]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq3 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 3]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq4 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 4]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq5 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 5]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               })
      ) 
    
    y %>%
      inner_join(., reference, "phon")
  }) %>%
  mutate(seq2 = sapply(seq2, unique))

chi_tokens_converted <- chi_na_id_sec %>%
  filter(str_detect(phon, "_")) %>%
  mutate(phon = phon %>%
           str_replace_all("UU", "AH")) %>%
  (function(y) {
    reference <- y[!duplicated(y$phon),] %>%
      mutate(split = str_split(phon, "_") %>%
               sapply(function(x) {
                 adj_comb_complete(x)
               })) %>%
      mutate(seq2 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 2]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq3 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 3]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq4 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 4]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq5 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 5]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               })
      ) 
    
    y %>%
      inner_join(., reference, "phon")
  }) %>%
  select(id.x:phon, split:seq5) %>%
  (function(x) {
    colnames(x) <- str_remove(colnames(x), "[.]{1}.*$")
    x
  }) %>%
  mutate(seq2 = sapply(seq2, unique)) %>%
  select(phon:seq5)

#### local funs ####
add_frame_n_words <- function(df, seq_num, name_var1, name_var2, temp_df = "mc") {
  if (temp_df == "cpwd") {
    temp <- cpwd_tokens %>%
      select(phon, seq_num) %>%
      (function(x) {
        x[[seq_num]][lengths(x[[seq_num]]) == 0] <- NA
        x
      }) %>%
      .[!duplicated(.$phon),] %>%
      apply(1, function(x) {
        tibble(phon = x$phon, seq = unlist(x[[seq_num]]))
      }) %>% 
      rbindlist() %>%
      na.omit() %>%
      arrange(phon)
  } else {
    temp <- mot_tokens_converted %>%
      select(phon, seq_num) %>%
      rbind(chi_tokens_converted %>%
              select(phon, seq_num)) %>%
      (function(x) {
        x[[seq_num]][lengths(x[[seq_num]]) == 0] <- NA
        x
      }) %>%
      .[!duplicated(.$phon),] %>%
      apply(1, function(x) {
        tibble(phon = x$phon, seq = unlist(x[[seq_num]]))
      }) %>% 
      rbindlist() %>%
      na.omit() %>%
      arrange(phon)
  }
  
  df %>%
    mutate(!!name_var1 := sapply(1:nrow(df), function(i) {
      filtered_temp <- temp %>%
        filter(phon != df$phon[i])
      
      temp1 <- str_count(filtered_temp %>%
                           .$seq,
                         paste("^", df[[seq_num]][[i]], "$", collapse = "|", sep="")) 
      
      temp1 %>% sum()
    }),
    !!name_var2 := sapply(1:nrow(df), function(i) {
      filtered_temp <- temp %>%
        filter(phon != df$phon[i])
      
      temp1 <- str_count(filtered_temp %>%
                           .$seq,
                         paste("^", df[[seq_num]][[i]], "$", collapse = "|", sep="")) 
      
      temp1 %>%
      {filtered_temp$phon[which(. == 1)] %>%
          unique() %>% sort()} %>%
        paste(collapse = " | ")
    }))
}

add_frame_freq <- function(df, seq_num, name_var1, temp_df = "mc") {
  if (temp_df == "cpwd") {
    temp <- cpwd_tokens %>%
      select(phon, seq_num) %>%
      (function(x) {
        x[[seq_num]][lengths(x[[seq_num]]) == 0] <- NA
        x
      }) %>%
      apply(1, function(x) {
        tibble(phon = x$phon, seq = unlist(x[[seq_num]]))
      }) %>% 
      rbindlist() %>%
      na.omit() %>%
      arrange(phon)
  } else {
    temp <- mot_tokens_converted %>%
      select(phon, seq_num) %>%
      rbind(chi_tokens_converted %>%
              select(phon, seq_num)) %>%
      (function(x) {
        x[[seq_num]][lengths(x[[seq_num]]) == 0] <- NA
        x
      }) %>%
      apply(1, function(x) {
        tibble(phon = x$phon, seq = unlist(x[[seq_num]]))
      }) %>% 
      rbindlist() %>%
      na.omit() %>%
      arrange(phon)
  }
  
  df %>%
    mutate(!!name_var1 := sapply(1:nrow(df), function(i) {
      filtered_temp <- temp %>%
        filter(phon != df$phon[i])
      
      temp1 <- str_count(filtered_temp %>%
                           .$seq,
                         paste("^", df[[seq_num]][[i]], "$", collapse = "|", sep="")) 
      
      temp1 %>% sum()
    }))
}

#### tables of frames vars ####
# manchester corpus nouns
mc_nouns <- extract_gram_types(chi_uni_cat, "N") %>%
  rbind(mot_nouns) %>%
  filter(!str_detect(word, "@")) %>%
  filter(str_detect(phon, "_")) %>%
  mutate(phon = phon %>%
           str_replace_all("UU", "AH")) %>%
  seq2_5() %>%
  .[!duplicated(.$phon),] %>%
  arrange(phon) %>%
  select(-id, -section, -cat, -word) %>%
  mutate(seq2 = sapply(seq2, unique))

mc_nouns %<>%
  add_frame_n_words("seq3", "3_frames", "3_frame_words") %>%
  add_frame_n_words("seq4", "4_frames", "4_frame_words") %>%
  add_frame_n_words("seq5", "5_frames", "5_frame_words")

mc_nouns %<>%
  add_frame_freq("seq3", "tot_3_frame_freqs") %>%
  add_frame_freq("seq4", "tot_4_frame_freqs") %>%
  add_frame_freq("seq5", "tot_5_frame_freqs")

# import CPWD
cpwd_types <- "cpwd_freqs_phon.txt" %>%
  read_tsv(col_names = FALSE) %>%
  .[-1,] %>%
  .[-nrow(.),] %>%
  separate(X1, c("phon", "freq"), sep = "\\) ") %>%
  mutate(freq = freq %>%
           str_remove("\\)$") %>%
           as.numeric(),
         phon = phon %>%
           str_remove_all("[0-9]+|^\\(\\(") %>%
           str_replace_all(" ", "_")) %>%
  arrange(phon)

cpwd_tokens <- lapply(1:nrow(cpwd_types), function(i) {
  tibble(phon = rep(cpwd_types$phon[i], cpwd_types$freq[i]))
}) %>%
  rbindlist()

cpwd_types %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5() %>%
  mutate(seq2 = sapply(seq2, unique))

cpwd_tokens %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5() %>%
  mutate(seq2 = sapply(seq2, unique))

cpwd_types %<>%
  add_frame_n_words("seq3", "3_frames", "3_frame_words", temp_df = "cpwd") %>%
  add_frame_n_words("seq4", "4_frames", "4_frame_words", temp_df = "cpwd") %>%
  add_frame_n_words("seq5", "5_frames", "5_frame_words", temp_df = "cpwd")

cpwd_types %<>%
  add_frame_freq("seq3", "tot_3_frame_freqs", temp_df = "cpwd") %>%
  add_frame_freq("seq4", "tot_4_frame_freqs", temp_df = "cpwd") %>%
  add_frame_freq("seq5", "tot_5_frame_freqs", temp_df = "cpwd")

# joined df
joined_nouns <- mc_nouns$phon %>%
  intersect(cpwd_types %>%
              .$phon %>%
              unique()) %>%
  (function(x) {
    mc_nouns %>%
      filter(phon %in% x) %>%
      arrange(phon) %>%
      select(phon, `3_frames`:`tot_5_frame_freqs`) %>%
      (function(y) {
        colnames(y) <- str_replace(colnames(y), "$", "_mc")
        y
      }) %>%
      cbind(cpwd_types %>%
              filter(phon %in% x) %>%
              .[!duplicated(.$phon), ] %>%
              arrange(phon) %>%
              select(phon, `3_frames`:`tot_5_frame_freqs`) %>%
              (function(y) {
                colnames(y) <- str_replace(colnames(y), "$", "_cpwd")
                y
              })
      )
  }) %>%
  select(-phon_cpwd) %>%
  rename(phon = phon_mc) %>%
  mutate(temp = `3_frames_mc` + `3_frames_cpwd`) %>%
  arrange(desc(temp)) %>%
  select(-temp)

write_csv(mc_nouns %>%
            select(-(syllable:seq5)) %>%
            arrange(desc(`3_frames`)),
          "mc_nouns.csv")

write_csv(cpwd_types %>%
            select(-(freq:seq5)) %>%
            arrange(desc(`3_frames`)),
          "cpwd_types.csv")

write_csv(joined_nouns,
          "mc&cpwd_nouns.csv")

#### summary table MC ####
mc_cpwd <- read_csv("mc&cpwd_nouns.csv")

mc_words <- (mot_tokens_converted$phon %>%
               c(chi_tokens_converted$phon))

prova_mc_words <- tibble(words = mc_words)

summary_mc <- mc_cpwd %>%
  select(target = phon) %>%
  mutate(len_target = str_split(target, "_") %>% 
           sapply(function(x) {
             x %>% 
               length()
           }) %>%
           round(2)) %>%
  filter(len_target >= 3)

summary_mc %<>%
  mutate(phon = target) %>%
  seq2_5() %>% 
  select(target:len_target, seq3:seq5) %>%
  mutate(seq3 = sapply(seq3, unique)) %>%
  mutate(seq4 = sapply(seq4, unique)) %>%
  mutate(seq5 = sapply(seq5, unique)) %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$seq[[i]] <- c(x$seq3[[i]], x$seq4[[i]], x$seq5[[i]])
    }
    x
  }) %>%
  select(-(seq3:seq5)) %>%
  mutate(seq = sapply(seq, function(x) {
    if (length(x) == 0) {
      list(NA)
    } else {
      x
    }
  })) %>%
  apply(1, function(x) {
    tibble(target = x$target,
           len_target = x$len_target,
           seq = x$seq)
  }) %>%
  rbindlist() 

summary_mc %<>%
  mutate(seq = unlist(seq)) %>%
  mutate(phon = target) %>%
  inner_join(., mc_cpwd, by = "phon") %>%
  select(target:seq, `3_frame_words_mc`, `4_frame_words_mc`, `5_frame_words_mc`) %>%
  unite(other, 
        `3_frame_words_mc`, 
        `4_frame_words_mc`, 
        `5_frame_words_mc`, sep = " | ") %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$shared[i] <- str_detect(x$other[i], x$seq[i]) 
    }
    x
  }) %>%
  filter(shared == TRUE) %>%
  select(target, len_target) %>%
  .[!duplicated(.$target),] 

summary_mc %<>%
  mutate(n_target = n()) %>%
  select(target, n_target, len_target) %>%
  mutate(avg_tot_len_target = mean(len_target, na.rm = TRUE) %>% round(2),
         sd_tot_len_target = sd(len_target, na.rm = TRUE) %>% round(2),
         freq_target = target %>%
           sapply(function(x) {
             fmatch(mc_words, x) %>% sum(na.rm=TRUE)
           })) %>%
  mutate(freq_target = round((freq_target*1000000) / length(mc_words), 2)) %>%
  mutate(avg_tot_freq_target = mean(freq_target, na.rm = TRUE) %>% round(2),
         sd_tot_freq_target = sd(freq_target, na.rm = TRUE) %>% round(2)) 

summary_mc %<>%
  mutate(phon = target) %>%
  seq2_5() %>% 
  select(target:sd_tot_freq_target, seq3:seq5) %>%
  mutate(seq3 = sapply(seq3, unique)) %>%
  mutate(seq4 = sapply(seq4, unique)) %>%
  mutate(seq5 = sapply(seq5, unique)) %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$seq[[i]] <- c(x$seq3[[i]], x$seq4[[i]], x$seq5[[i]])
    }
    x
  }) %>%
  select(-(seq3:seq5)) %>%
  mutate(seq = sapply(seq, function(x) {
    if (length(x) == 0) {
      list(NA)
    } else {
      x
    }
  })) %>%
  apply(1, function(x) {
    tibble(target = x$target,
           n_target = x$n_target,
           len_target = x$len_target,
           avg_tot_len_target = x$avg_tot_len_target,
           sd_tot_len_target = x$sd_tot_len_target,
           freq_target = x$freq_target,
           avg_tot_freq_target = x$avg_tot_freq_target,
           sd_tot_freq_target = x$sd_tot_freq_target,
           seq = x$seq)
  }) %>%
  rbindlist() 

summary_mc %<>%
  mutate(seq = unlist(seq)) %>%
  mutate(phon = target) %>%
  inner_join(., mc_cpwd, by = "phon") %>%
  select(target:seq, `3_frame_words_mc`, `4_frame_words_mc`, `5_frame_words_mc`) %>%
  unite(other, 
        `3_frame_words_mc`, 
        `4_frame_words_mc`, 
        `5_frame_words_mc`, sep = " | ") %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$shared[i] <- str_detect(x$other[i], x$seq[i]) 
    }
    x
  }) %>%
  filter(shared == TRUE) %>%
  select(-shared)

summary_mc %<>%
  rename(frame = seq) %>%
  select(other, target:frame) %>%
  group_by(target) %>%
  mutate(n_frame = c(n())) %>%
  ungroup() 

summary_mc %<>%
  (function(x) {
    reference <- x %>%
      select(target, n_frame) %>%
      group_by(target) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      mutate(avg_tot_n_frame = round(mean(n_frame, na.rm = TRUE), 2),
             sd_tot_n_frame = round(sd(n_frame, na.rm = TRUE), 2))
    
    x %>%
      inner_join(., reference, by = c("target", "n_frame"))
  }) 

summary_mc %<>%
  mutate(len_frame = str_split(frame, "_") %>% 
           sapply(function(x) {
             x %>% 
               length()
           })) %>%
  group_by(target) %>%
  mutate(avg_len_frame = round(mean(len_frame, na.rm = TRUE), 2),
         sd_len_frame = round(sd(len_frame, na.rm = TRUE), 2)) %>%
  ungroup() %>%
  mutate(other = str_split(other, " \\| ") %>%
           lapply(unique)) %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$other[[i]] <- x$other[[i]][str_detect(x$other[[i]], x$frame[i])] %>% sort()
    } 
    x
  }) 

summary_mc %<>%
  (function(x) {
    reference <- x %>%
      select(target, avg_len_frame) %>%
      group_by(target) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      mutate(avg_tot_len_frame = round(mean(avg_len_frame, na.rm = TRUE), 2),
             sd_tot_len_frame = round(sd(avg_len_frame, na.rm = TRUE), 2))
    
    x %>%
      inner_join(., reference, by = c("target", "avg_len_frame"))
  })

# freq_frame 
summary_mc %<>%
  rowwise() %>%
  mutate(freq_frame = mc_words[which(mc_words %in% other)] %>% 
           length() %>%
           {((.*1000000) / length(mc_words)) %>%
               round(2)}) %>%
  ungroup() %>%
  group_by(target) %>%
  mutate(avg_freq_frame = mean(freq_frame) %>% round(2),
         sd_freq_frame = sd(freq_frame) %>% round(2)) %>%
  ungroup() 

summary_mc %<>%
  (function(x) {
    reference <- x %>%
      select(target, avg_freq_frame) %>%
      group_by(target) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      mutate(avg_tot_freq_frame = round(mean(avg_freq_frame, na.rm = TRUE), 2),
             sd_tot_freq_frame = round(sd(avg_freq_frame, na.rm = TRUE), 2))
    
    x %>%
      inner_join(., reference, by = c("target", "avg_freq_frame"))
  })

summary_mc %<>%
  apply(1, function(x) {
    tibble(target = x$target,
           n_target = x$n_target,
           len_target = x$len_target,
           avg_tot_len_target = x$avg_tot_len_target,
           sd_tot_len_target = x$sd_tot_len_target,
           freq_target = x$freq_target,
           avg_tot_freq_target = x$avg_tot_freq_target,
           sd_tot_freq_target = x$sd_tot_freq_target,
           frame = x$frame,
           n_frame = x$n_frame,
           avg_tot_n_frame = x$avg_tot_n_frame,
           sd_tot_n_frame = x$sd_tot_n_frame,
           len_frame = x$len_frame,
           avg_len_frame = x$avg_len_frame,
           sd_len_frame = x$sd_len_frame,
           avg_tot_len_frame = x$avg_tot_len_frame,
           sd_tot_len_frame = x$sd_tot_len_frame,
           freq_frame = x$freq_frame,
           avg_freq_frame = x$avg_freq_frame,
           sd_freq_frame = x$sd_freq_frame,
           avg_tot_freq_frame = x$avg_tot_freq_frame,
           sd_tot_freq_frame = x$sd_tot_freq_frame,
           other = unlist(x$other))
  }) %>%
  rbindlist()

# freq_other
summary_mc %<>%
  group_by(target, frame) %>%
  mutate(n_other = n()) %>%
  ungroup() 

summary_mc %<>%
  (function(x) {
    reference <- summary_mc %>%
      select(target, frame, n_other) %>%
      group_by(target, frame) %>%
      filter(row_number() == 1) %>%
      ungroup() %>%
      group_by(target) %>%
      summarise(avg_n_other = mean(n_other, na.rm = TRUE) %>% round(2),
                sd_n_other = sd(n_other, na.rm = TRUE) %>% round(2)) %>%
      ungroup()
    
    x %>%
      inner_join(., reference, by = "target")
  })

summary_mc %<>%
  (function(x) {
    reference <- x %>%
      select(target, avg_n_other) %>%
      group_by(target) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      mutate(avg_tot_n_other = round(mean(avg_n_other, na.rm = TRUE), 2),
             sd_tot_n_other = round(sd(avg_n_other, na.rm = TRUE), 2))
    
    x %>%
      inner_join(., reference, by = c("target", "avg_n_other"))
  })

summary_mc %<>%
  (function(x) {
    reference <- x %>%
      select(other) %>%
      distinct(other) %>%
      mutate(freq_other = other %>%
               sapply(function(y) {
                 fmatch(mc_words, y) %>% sum(na.rm=TRUE)
               }))
    
    x %>%
      inner_join(., reference, by = "other")
  }) %>%
  mutate(freq_other = ((freq_other*1000000) / length(mc_words)) %>%
           round(2)) 

summary_mc %<>%
  (function(x) {
    reference <- x %>%
      select(target, frame, freq_other) %>%
      group_by(target, frame) %>%
      summarise(m_freq_other = mean(freq_other, na.rm = TRUE) %>% round(2)) %>%
      ungroup() %>%
      group_by(target) %>%
      summarise(avg_freq_other = mean(m_freq_other, na.rm = TRUE) %>% round(2),
                sd_freq_other = sd(m_freq_other, na.rm = TRUE) %>% round(2))
    
    x %>%
      inner_join(., reference, by = "target")
  })

summary_mc %<>%
  (function(x) {
    reference <- x %>%
      select(target, avg_freq_other) %>%
      group_by(target, avg_freq_other) %>%
      filter(row_number() == 1) %>%
      ungroup() %>%
      mutate(avg_tot_freq_other = mean(avg_freq_other, na.rm = TRUE) %>% round(2),
             sd_tot_freq_other = sd(avg_freq_other, na.rm = TRUE) %>% round(2)) 
    
    x %>%
      inner_join(., reference, by = c("target", "avg_freq_other"))
  })

summary_mc %<>%
  select(target, freq_target, frame, freq_frame, other, freq_other,
         n_target, n_frame, n_other,
         len_target, len_frame,
         avg_tot_freq_target, sd_tot_freq_target,
         avg_freq_frame, sd_freq_frame,
         avg_tot_freq_frame, sd_tot_freq_frame,
         avg_freq_other, sd_freq_other,
         avg_tot_freq_other, sd_tot_freq_other,
         avg_tot_n_frame, sd_tot_n_frame,
         avg_n_other, sd_n_other,
         avg_tot_n_other, sd_tot_n_other,
         avg_tot_len_target, sd_tot_len_target,
         avg_len_frame, sd_len_frame,
         avg_tot_len_frame, sd_tot_len_frame)

summary_mc_info <- summary_mc %>%
  select(contains("_tot_")) %>%
  .[1,] %>%
  (function(x) {
    tibble(var = colnames(x) %>% str_remove("^sd_tot_|^avg_tot_") %>% unique(),
           avg_tot = x %>% 
             select(contains("avg_")) %>%
             unlist(),
           sd_tot = x %>% 
             select(contains("sd_")) %>%
             unlist())
  })

summary_mc %<>%
  select(-contains("_tot_"))

write_csv(summary_mc, "summary_mc.csv")
write_csv(summary_mc_info, "summary_mc_info.csv")
rm(prova_mc_words, mc_words)


#### summary_table CPWD ####
cpwd_words <- cpwd_tokens$phon

prova_cpwd_words <- tibble(words = cpwd_words)

summary_cpwd <- mc_cpwd %>%
  select(target = phon) %>%
  mutate(len_target = str_split(target, "_") %>% 
           sapply(function(x) {
             x %>% 
               length()
           }) %>%
           round(2)) %>%
  filter(len_target >= 3)

summary_cpwd %<>%
  mutate(phon = target) %>%
  seq2_5() %>% 
  select(target:len_target, seq3:seq5) %>%
  mutate(seq3 = sapply(seq3, unique)) %>%
  mutate(seq4 = sapply(seq4, unique)) %>%
  mutate(seq5 = sapply(seq5, unique)) %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$seq[[i]] <- c(x$seq3[[i]], x$seq4[[i]], x$seq5[[i]])
    }
    x
  }) %>%
  select(-(seq3:seq5)) %>%
  mutate(seq = sapply(seq, function(x) {
    if (length(x) == 0) {
      list(NA)
    } else {
      x
    }
  })) %>%
  apply(1, function(x) {
    tibble(target = x$target,
           len_target = x$len_target,
           seq = x$seq)
  }) %>%
  rbindlist() 

summary_cpwd %<>%
  mutate(seq = unlist(seq)) %>%
  mutate(phon = target) %>%
  inner_join(., mc_cpwd, by = "phon") %>%
  select(target:seq, `3_frame_words_cpwd`, `4_frame_words_cpwd`, `5_frame_words_cpwd`) %>%
  unite(other, 
        `3_frame_words_cpwd`, 
        `4_frame_words_cpwd`, 
        `5_frame_words_cpwd`, sep = " | ") %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$shared[i] <- str_detect(x$other[i], x$seq[i]) 
    }
    x
  }) %>%
  filter(shared == TRUE) %>%
  select(target, len_target) %>%
  .[!duplicated(.$target),] 

summary_cpwd %<>%
  mutate(n_target = n()) %>%
  select(target, n_target, len_target) %>%
  mutate(avg_tot_len_target = mean(len_target, na.rm = TRUE) %>% round(2),
         sd_tot_len_target = sd(len_target, na.rm = TRUE) %>% round(2),
         freq_target = target %>%
           sapply(function(x) {
             fmatch(cpwd_words, x) %>% sum(na.rm=TRUE)
           })) %>%
  #mutate(freq_target = round((freq_target*1000000) / length(cpwd_words), 2)) %>%
  mutate(avg_tot_freq_target = mean(freq_target, na.rm = TRUE) %>% round(2),
         sd_tot_freq_target = sd(freq_target, na.rm = TRUE) %>% round(2)) 

summary_cpwd %<>%
  mutate(phon = target) %>%
  seq2_5() %>% 
  select(target:sd_tot_freq_target, seq3:seq5) %>%
  mutate(seq3 = sapply(seq3, unique)) %>%
  mutate(seq4 = sapply(seq4, unique)) %>%
  mutate(seq5 = sapply(seq5, unique)) %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$seq[[i]] <- c(x$seq3[[i]], x$seq4[[i]], x$seq5[[i]])
    }
    x
  }) %>%
  select(-(seq3:seq5)) %>%
  mutate(seq = sapply(seq, function(x) {
    if (length(x) == 0) {
      list(NA)
    } else {
      x
    }
  })) %>%
  apply(1, function(x) {
    tibble(target = x$target,
           n_target = x$n_target,
           len_target = x$len_target,
           avg_tot_len_target = x$avg_tot_len_target,
           sd_tot_len_target = x$sd_tot_len_target,
           freq_target = x$freq_target,
           avg_tot_freq_target = x$avg_tot_freq_target,
           sd_tot_freq_target = x$sd_tot_freq_target,
           seq = x$seq)
  }) %>%
  rbindlist() 

summary_cpwd %<>%
  mutate(seq = unlist(seq)) %>%
  mutate(phon = target) %>%
  inner_join(., mc_cpwd, by = "phon") %>%
  select(target:seq, `3_frame_words_cpwd`, `4_frame_words_cpwd`, `5_frame_words_cpwd`) %>%
  unite(other, 
        `3_frame_words_cpwd`, 
        `4_frame_words_cpwd`, 
        `5_frame_words_cpwd`, sep = " | ") %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$shared[i] <- str_detect(x$other[i], x$seq[i]) 
    }
    x
  }) %>%
  filter(shared == TRUE) %>%
  select(-shared)

summary_cpwd %<>%
  rename(frame = seq) %>%
  select(other, target:frame) %>%
  group_by(target) %>%
  mutate(n_frame = c(n())) %>%
  ungroup() 

summary_cpwd %<>%
  (function(x) {
    reference <- x %>%
      select(target, n_frame) %>%
      group_by(target) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      mutate(avg_tot_n_frame = round(mean(n_frame, na.rm = TRUE), 2),
             sd_tot_n_frame = round(sd(n_frame, na.rm = TRUE), 2))
    
    x %>%
      inner_join(., reference, by = c("target", "n_frame"))
  }) 

summary_cpwd %<>%
  mutate(len_frame = str_split(frame, "_") %>% 
           sapply(function(x) {
             x %>% 
               length()
           })) %>%
  group_by(target) %>%
  mutate(avg_len_frame = round(mean(len_frame, na.rm = TRUE), 2),
         sd_len_frame = round(sd(len_frame, na.rm = TRUE), 2)) %>%
  ungroup() %>%
  mutate(other = str_split(other, " \\| ") %>%
           lapply(unique)) %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$other[[i]] <- x$other[[i]][str_detect(x$other[[i]], x$frame[i])] %>% sort()
    } 
    x
  }) 

summary_cpwd %<>%
  (function(x) {
    reference <- x %>%
      select(target, avg_len_frame) %>%
      group_by(target) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      mutate(avg_tot_len_frame = round(mean(avg_len_frame, na.rm = TRUE), 2),
             sd_tot_len_frame = round(sd(avg_len_frame, na.rm = TRUE), 2))
    
    x %>%
      inner_join(., reference, by = c("target", "avg_len_frame"))
  })

# freq_frame 
summary_cpwd %<>%
  rowwise() %>%
  mutate(freq_frame = cpwd_words[which(cpwd_words %in% other)] %>% 
           length()) %>%
  ungroup() %>%
  group_by(target) %>%
  mutate(avg_freq_frame = mean(freq_frame) %>% round(2),
         sd_freq_frame = sd(freq_frame) %>% round(2)) %>%
  ungroup() 

summary_cpwd %<>%
  (function(x) {
    reference <- x %>%
      select(target, avg_freq_frame) %>%
      group_by(target) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      mutate(avg_tot_freq_frame = round(mean(avg_freq_frame, na.rm = TRUE), 2),
             sd_tot_freq_frame = round(sd(avg_freq_frame, na.rm = TRUE), 2))
    
    x %>%
      inner_join(., reference, by = c("target", "avg_freq_frame"))
  })

summary_cpwd %<>%
  apply(1, function(x) {
    tibble(target = x$target,
           n_target = x$n_target,
           len_target = x$len_target,
           avg_tot_len_target = x$avg_tot_len_target,
           sd_tot_len_target = x$sd_tot_len_target,
           freq_target = x$freq_target,
           avg_tot_freq_target = x$avg_tot_freq_target,
           sd_tot_freq_target = x$sd_tot_freq_target,
           frame = x$frame,
           n_frame = x$n_frame,
           avg_tot_n_frame = x$avg_tot_n_frame,
           sd_tot_n_frame = x$sd_tot_n_frame,
           len_frame = x$len_frame,
           avg_len_frame = x$avg_len_frame,
           sd_len_frame = x$sd_len_frame,
           avg_tot_len_frame = x$avg_tot_len_frame,
           sd_tot_len_frame = x$sd_tot_len_frame,
           freq_frame = x$freq_frame,
           avg_freq_frame = x$avg_freq_frame,
           sd_freq_frame = x$sd_freq_frame,
           avg_tot_freq_frame = x$avg_tot_freq_frame,
           sd_tot_freq_frame = x$sd_tot_freq_frame,
           other = unlist(x$other))
  }) %>%
  rbindlist()

# freq_other
summary_cpwd %<>%
  group_by(target, frame) %>%
  mutate(n_other = n()) %>%
  ungroup() 

summary_cpwd %<>%
  (function(x) {
    reference <- summary_cpwd %>%
      select(target, frame, n_other) %>%
      group_by(target, frame) %>%
      filter(row_number() == 1) %>%
      ungroup() %>%
      group_by(target) %>%
      summarise(avg_n_other = mean(n_other, na.rm = TRUE) %>% round(2),
                sd_n_other = sd(n_other, na.rm = TRUE) %>% round(2)) %>%
      ungroup()
    
    x %>%
      inner_join(., reference, by = "target")
  })

summary_cpwd %<>%
  (function(x) {
    reference <- x %>%
      select(target, avg_n_other) %>%
      group_by(target) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      mutate(avg_tot_n_other = round(mean(avg_n_other, na.rm = TRUE), 2),
             sd_tot_n_other = round(sd(avg_n_other, na.rm = TRUE), 2))
    
    x %>%
      inner_join(., reference, by = c("target", "avg_n_other"))
  })

summary_cpwd %<>%
  (function(x) {
    reference <- x %>%
      select(other) %>%
      distinct(other) %>%
      mutate(freq_other = other %>%
               sapply(function(y) {
                 fmatch(cpwd_words, y) %>% sum(na.rm=TRUE)
               }))
    
    x %>%
      inner_join(., reference, by = "other")
  }) 

summary_cpwd %<>%
  (function(x) {
    reference <- x %>%
      select(target, frame, freq_other) %>%
      group_by(target, frame) %>%
      summarise(m_freq_other = mean(freq_other, na.rm = TRUE) %>% round(2)) %>%
      ungroup() %>%
      group_by(target) %>%
      summarise(avg_freq_other = mean(m_freq_other, na.rm = TRUE) %>% round(2),
                sd_freq_other = sd(m_freq_other, na.rm = TRUE) %>% round(2))
    
    x %>%
      inner_join(., reference, by = "target")
  })

summary_cpwd %<>%
  (function(x) {
    reference <- x %>%
      select(target, avg_freq_other) %>%
      group_by(target, avg_freq_other) %>%
      filter(row_number() == 1) %>%
      ungroup() %>%
      mutate(avg_tot_freq_other = mean(avg_freq_other, na.rm = TRUE) %>% round(2),
             sd_tot_freq_other = sd(avg_freq_other, na.rm = TRUE) %>% round(2)) 
    
    x %>%
      inner_join(., reference, by = c("target", "avg_freq_other"))
  })

summary_cpwd %<>%
  select(target, freq_target, frame, freq_frame, other, freq_other,
         n_target, n_frame, n_other,
         len_target, len_frame,
         avg_tot_freq_target, sd_tot_freq_target,
         avg_freq_frame, sd_freq_frame,
         avg_tot_freq_frame, sd_tot_freq_frame,
         avg_freq_other, sd_freq_other,
         avg_tot_freq_other, sd_tot_freq_other,
         avg_tot_n_frame, sd_tot_n_frame,
         avg_n_other, sd_n_other,
         avg_tot_n_other, sd_tot_n_other,
         avg_tot_len_target, sd_tot_len_target,
         avg_len_frame, sd_len_frame,
         avg_tot_len_frame, sd_tot_len_frame)

summary_cpwd_info <- summary_cpwd %>%
  select(contains("_tot_")) %>%
  .[1,] %>%
  (function(x) {
    tibble(var = colnames(x) %>% str_remove("^sd_tot_|^avg_tot_") %>% unique(),
           avg_tot = x %>% 
             select(contains("avg_")) %>%
             unlist(),
           sd_tot = x %>% 
             select(contains("sd_")) %>%
             unlist())
  })

summary_cpwd %<>%
  select(-contains("_tot_"))

write_csv(summary_cpwd, "summary_cpwd.csv")
write_csv(summary_cpwd_info, "summary_cpwd_info.csv")
rm(prova_cpwd_words, cpwd_words)






summary_mc_info %>% 
  select(contains("avg_")) %>%
  unlist() %>%
  class()
