# load libraries
lib <- c("magrittr", "tidyverse", 
         "data.table", "fastmatch",
         "beepr", "mailR", "janitor")
lapply(lib, require, character.only = TRUE)
rm(lib)

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
  select(-id, -section, -cat, -word) 

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
  seq2_5() 

cpwd_tokens %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5() 

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