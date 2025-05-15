run_simulation_st <- function(country_analysis, n_i, p_i_mild, p_i_moderate, r_i_MV, c_trt, c_MV, c_ICU, c_GW, c_REH, c_REC, c_D, hr_mild_moderate, hr_moderate_severe, hr_mild_severe, hr_MV_ICU_mild, hr_MV_ICU_moderate, hr_MV_ICU_severe){
  
  input_parameters <- function(country_analysis, n_i, p_i_mild, p_i_moderate, r_i_MV, c_trt, c_MV, c_ICU, c_GW, c_REH, c_REC, c_D, hr_mild_moderate, hr_moderate_severe, hr_mild_severe, hr_MV_ICU_mild, hr_MV_ICU_moderate, hr_MV_ICU_severe){
    country_analysis <- country_analysis
    cycle_length   <- 1/365  # cycle length equal to one year
    n_cycles       <- 60    # time horizon, number of cycles
    v_names_states <- c("MV_mild",
                        "ICU_mild",       
                        "GW_mild",        
                        "MV_moderate",
                        "ICU_moderate",   
                        "GW_moderate",    
                        "MV_severe",
                        "ICU_severe",      
                        "GW_severe",
                        "REH",
                        "REC",
                        "D")
    
    n_states       <- length(v_names_states) # number of health states 
    v_names_str   <- c("Standard of care",   # store the strategy names
                       "FX06 treatment") 
    n_str         <- length(v_names_str)     # number of strategies
    
    # Load functions and pre-calculated data
    df_p_mv_icu <- read.csv("https://raw.githubusercontent.com/clazinusveijer/shepards_bia/refs/heads/main/model_input_data/df_p_mv_icu.csv", sep = ',')
    
    #load the file with hazard ratios for return to work per ICU length of stay (no. days)
    df_hr_iculos_rtw <- read.csv("https://raw.githubusercontent.com/clazinusveijer/shepards_bia/refs/heads/main/model_input_data/df_hr_iculos_rtw.csv", sep = ',')
    #extend the file with hazard ratio of the final day until to total number of cycles
    df_hr_iculos_rtw <- rbind(df_hr_iculos_rtw, data.frame(n_cycles_ICU = rep(max(df_hr_iculos_rtw$n_cycles_ICU+1):n_cycles), hr_RTW = df_hr_iculos_rtw$hr_RTW[which(df_hr_iculos_rtw$n_cycles_ICU == 30)]))
    
    p_mort <- read.csv("https://raw.githubusercontent.com/clazinusveijer/shepards_bia/refs/heads/main/model_input_data/df_p_mort_cy.csv", sep = ',')
    df_p_REC_D <- p_mort %>% filter(country == country_analysis)
    
    n_i            <- n_i 
    p_i_mild       <- p_i_mild 
    p_i_moderate   <- p_i_moderate 
    #p_i_severe     <- p_i_severe
    n_i_mild       <- round(n_i * p_i_mild,0)
    n_i_moderate   <- round(n_i * p_i_moderate,0)
    n_i_severe     <- n_i - n_i_mild - n_i_moderate #round(n_i * p_i_severe,0)
    
    # initial number of patients per health state
    r_i_MV  <- r_i_MV #0.856
    r_i_ICU <- 1-r_i_MV # 0.144 
    n_i_mild_MV  <- round(n_i_mild*r_i_MV,0)
    n_i_mild_ICU <- n_i_mild - n_i_mild_MV #round(n_i_mild*r_i_ICU,0)
    n_i_moderate_MV  <- round(n_i_moderate*r_i_MV,0)
    n_i_moderate_ICU <- n_i_moderate - n_i_moderate_MV #round(n_i_moderate*r_i_ICU,0)
    n_i_severe_MV  <- round(n_i_severe*r_i_MV,0)
    n_i_severe_ICU <- n_i_severe - n_i_severe_MV #round(n_i_severe*r_i_ICU,0)
    
    hr_mild_moderate        <- hr_mild_moderate  
    hr_moderate_severe      <- hr_moderate_severe
    hr_mild_severe          <- hr_mild_severe
    
    hr_MV_ICU_mild          <- hr_MV_ICU_mild
    hr_MV_ICU_moderate      <- hr_MV_ICU_moderate
    hr_MV_ICU_severe        <- hr_MV_ICU_severe
    
    # Disease progression up to 7 days post-diagnosis of ARDS
    d_progression_period    <- 7
    r_mild_moderate         <- 0.258
    r_moderate_severe       <- 0.127
    r_mild_severe           <- 0.045
    p_mild_moderate_soc     <- 1-((1-r_mild_moderate)^(1/d_progression_period))
    p_moderate_severe_soc   <- 1-((1-r_moderate_severe)^(1/d_progression_period))
    p_mild_severe_soc       <- 1-((1-r_mild_severe)^(1/d_progression_period))
    
    p_mild_moderate_soc     <- c(p_mild_moderate_soc, 0)
    p_moderate_severe_soc   <- c(p_moderate_severe_soc, 0)
    p_mild_severe_soc       <- c(p_mild_severe_soc, 0)
    
    r_mild_moderate_trt     <- r_mild_moderate * hr_mild_moderate
    r_moderate_severe_trt   <- r_moderate_severe *  hr_moderate_severe
    r_mild_severe_trt       <- r_mild_severe * hr_mild_severe
    
    p_mild_moderate_trt     <- 1-((1-r_mild_moderate_trt)^(1/d_progression_period))
    p_moderate_severe_trt   <- 1-((1-r_moderate_severe_trt)^(1/d_progression_period))
    p_mild_severe_trt       <- 1-((1-r_mild_severe_trt)^(1/d_progression_period))
    
    p_mild_moderate_trt     <- c(p_mild_moderate_trt, 0)
    p_moderate_severe_trt   <- c(p_moderate_severe_trt, 0)
    p_mild_severe_trt       <- c(p_mild_severe_trt, 0)
    
    # 28-day mortality
    r_mort_28_mild          <- 0.296
    r_mort_28_moderate      <- 0.352
    r_mort_28_severe        <- 0.409
    p_D_mild                <- 1 - (1 - r_mort_28_mild)^(1 / 28)
    p_D_moderate            <- 1 - (1 - r_mort_28_moderate)^(1 / 28)
    p_D_severe              <- 1 - (1 - r_mort_28_severe)^(1 / 28)
    
    ### MV LOS
    median_MV_LOS_mild      <- 6
    median_MV_LOS_moderate  <- 8
    median_MV_LOS_severe    <- 11
    
    ### time-dependent MV LOS
    p_MV_ICU_mild_soc         <- df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'mild')]
    p_MV_ICU_moderate_soc     <- df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'moderate')]
    p_MV_ICU_severe_soc       <- df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'severe')]
    
    p_MV_ICU_mild_trt         <- df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'mild')]*hr_MV_ICU_mild
    p_MV_ICU_severe_trt       <- df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'severe')]*hr_MV_ICU_moderate
    p_MV_ICU_moderate_trt     <- df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'moderate')]*hr_MV_ICU_severe
    
    # ICU LOS
    median_ICU_LOS_mild     <- 10
    median_ICU_LOS_moderate <- 12
    median_ICU_LOS_severe   <- 14
    median_ICU_LOS_mild_net     <- median_ICU_LOS_mild - median_MV_LOS_mild
    median_ICU_LOS_moderate_net <- median_ICU_LOS_moderate - median_MV_LOS_moderate
    median_ICU_LOS_severe_net   <- median_ICU_LOS_severe - median_MV_LOS_severe
    
    r_ICU_mild              <- -log(0.5)/median_ICU_LOS_mild_net
    r_ICU_moderate          <- -log(0.5)/median_ICU_LOS_moderate_net
    r_ICU_severe            <- -log(0.5)/median_ICU_LOS_severe_net
    
    p_ICU_GW_mild           <- 1 - exp(-r_ICU_mild)
    p_ICU_GW_moderate       <- 1 - exp(-r_ICU_moderate)
    p_ICU_GW_severe         <- 1 - exp(-r_ICU_severe)
    
    # GW LOS
    median_GW_LOS_mild      <- 23
    median_GW_LOS_moderate  <- 22
    median_GW_LOS_severe    <- 26
    median_GW_LOS_mild_net      <- median_GW_LOS_mild - median_ICU_LOS_mild
    median_GW_LOS_moderate_net  <- median_GW_LOS_moderate - median_ICU_LOS_moderate
    median_GW_LOS_severe_net    <- median_GW_LOS_severe - median_ICU_LOS_severe
    
    r_GW_mild               <- -log(0.5)/median_GW_LOS_mild_net
    r_GW_moderate           <- -log(0.5)/median_GW_LOS_moderate_net
    r_GW_severe             <- -log(0.5)/median_GW_LOS_severe_net
    
    p_GW_REH_mild           <- 1 - exp(-r_GW_mild)
    p_GW_REH_moderate       <- 1 - exp(-r_GW_moderate)
    p_GW_REH_severe         <- 1 - exp(-r_GW_severe)
    
    LTS_RTW_median_d <- 203
    p_REH_REC     <- 1-(exp(-(-log(0.5)/LTS_RTW_median_d)))
    
    p_REH_D_st_1     <- 1-((1-0.042)^(1/(2*365))) #related to p_REH_REC and specific for age <65
    p_REH_D_st_2     <- 1-(1-0.356)^(1/(2*30)) #difference in mortality between day 30 and day 365, not related to p_REH_REC and specific for age >=65 
    df_p_REH_D <- data.frame(age_group_id = c(1,2), p_REH_D = c(p_REH_D_st_1, p_REH_D_st_2))
    
    #### Customized costs 
    c_MV      <- c_MV  # Each additional day on mechanically ventilated ICU patients
    c_ICU     <- c_ICU # Each additional day on non-mechanically ventilated ICU patients
    c_GW      <- c_GW  # Each additional day on the general ward for non-ICU patients
    c_REH     <- c_REH
    c_REC     <- c_REC
    c_D       <- c_D     # annual cost of being dead
    c_trt     <- c_trt # five-day treatment 600/5
    df_c      <- data.frame(type = c("c_MV", "c_ICU", "c_GW", "c_REH", "c_REC", "c_D", "c_trt"), value = c(c_MV, c_ICU, c_GW, c_REH, c_REC, c_D, c_trt))
    
    # sample from age distribution an initial age for every individual
    mean_age <- 61.5
    sd_age <- 14.9
    v_age_init <- rnorm(n_i, mean = mean_age, sd = sd_age)
    min_age <- 16
    max_age <- 99
    df_age <- data.frame(age = round(v_age_init, 0))
    df_age <- df_age %>% mutate(age = if_else(age<min_age, min_age, if_else(age>max_age, max_age, age)))
    df_age <- df_age %>% group_by(age) %>% summarize(prop = n()/n_i)
    v_age_init  <- sample(x = df_age$age, prob = df_age$prop, size = n_i, replace = TRUE) 
    
    v_M_init_mild <- rep(c("MV_mild", "ICU_mild"), times = c(n_i_mild_MV, n_i_mild_ICU))
    df_X_mild <- data.frame(ID = 1:n_i_mild, Severity = "mild", M_init = v_M_init_mild, n_cycles_MV = if_else(v_M_init_mild == "MV_mild", 1, 0)
                            ,   n_cycles_ICU = 1, trt_YN = 1, mild_moderate_YN = 1, mild_severe_YN = 1, moderate_severe_YN = 2)
    v_M_init_moderate <- rep(c("MV_moderate", "ICU_moderate"), times = c(n_i_moderate_MV, n_i_moderate_ICU))
    df_X_moderate <- data.frame(ID = (n_i_mild+1):(n_i_mild+n_i_moderate), Severity = "moderate", M_init = v_M_init_moderate, n_cycles_MV = if_else(v_M_init_moderate == "MV_moderate", 1, 0)
                                , n_cycles_ICU = 1, trt_YN = 1, mild_moderate_YN = 2, mild_severe_YN = 2, moderate_severe_YN =   1)
    v_M_init_severe <- rep(c("MV_severe", "ICU_severe"), times = c(n_i_severe_MV, n_i_severe_ICU))
    df_X_severe <- data.frame(ID = (n_i_mild+n_i_moderate+1):n_i, Severity = "severe", M_init = v_M_init_severe, n_cycles_MV = if_else(v_M_init_severe == "MV_severe", 1, 0)
                              , n_cycles_ICU = 1, trt_YN = 1, mild_moderate_YN = 2, mild_severe_YN = 2, moderate_severe_YN = 2)
    
    df_X <- rbind(df_X_mild, df_X_moderate, df_X_severe)
    df_X <- df_X %>% mutate(Age = v_age_init,
                            age_group_id = if_else(Age < 65, 1, 2)
    )
    
    inputs <- list(country_analysis = country_analysis
                   , cycle_length = cycle_length
                   , n_cycles = n_cycles
                   , v_names_states = v_names_states
                   , v_names_str = v_names_str
                   , n_str = n_str
                   , n_states = n_states
                   , n_i = n_i
                   , df_hr_iculos_rtw = df_hr_iculos_rtw
                   , df_p_REH_D = df_p_REH_D
                   , df_p_REC_D = df_p_REC_D 
                   , df_c = df_c
                   , df_X = df_X
                   , p_mild_moderate_soc = p_mild_moderate_soc
                   , p_moderate_severe_soc = p_moderate_severe_soc
                   , p_mild_severe_soc = p_mild_severe_soc
                   , p_mild_moderate_trt = p_mild_moderate_trt
                   , p_moderate_severe_trt = p_moderate_severe_trt
                   , p_mild_severe_trt = p_mild_severe_trt
                   , p_D_mild = p_D_mild
                   , p_D_moderate = p_D_moderate
                   , p_D_severe = p_D_severe
                   , p_MV_ICU_mild_soc = p_MV_ICU_mild_soc
                   , p_MV_ICU_moderate_soc = p_MV_ICU_moderate_soc
                   , p_MV_ICU_severe_soc = p_MV_ICU_severe_soc
                   , p_MV_ICU_mild_trt = p_MV_ICU_mild_trt
                   , p_MV_ICU_moderate_trt = p_MV_ICU_moderate_trt
                   , p_MV_ICU_severe_trt = p_MV_ICU_severe_trt
                   , p_ICU_GW_mild = p_ICU_GW_mild
                   , p_ICU_GW_moderate = p_ICU_GW_moderate
                   , p_ICU_GW_severe = p_ICU_GW_severe
                   , p_GW_REH_mild = p_GW_REH_mild
                   , p_GW_REH_moderate = p_GW_REH_moderate
                   , p_GW_REH_severe = p_GW_REH_severe
                   , p_REH_REC = p_REH_REC
    )
    return(inputs)
  }
  
  samplev <- function(m_Probs, m = 1) {
    lev <- dimnames(m_Probs)[[2]]  # extract the names of the health states considered for sampling
    n_samp <- nrow(m_Probs)
    u <- runif(n_samp, min = 0, max = 1)
    v_sum_p <- matrixStats::rowCumsums(m_Probs)
    v_cat <- lev[max.col(v_sum_p >= u, ties.method = "first")]
    return(v_cat)
  }
  
  Probs <- function(M_t, t, Trt, input_params, df_X, df_hr_iculos_rtw, df_p_REH_D, df_p_REC_D) { 
    #df_X <- as.data.frame(input_params$df_X)
    # Treatment specific transition probabilities
    if (Trt == "Standard of care") {
      p_mild_moderate <- input_params$p_mild_moderate_soc
    } else if (Trt == "FX06 treatment") {
      p_mild_moderate <- input_params$p_mild_moderate_trt
    }  else { 
      warning("Invalid treatment type (mild_moderate)") 
    }
    if (Trt == "Standard of care") {
      p_moderate_severe <- input_params$p_moderate_severe_soc
    } else if (Trt == "FX06 treatment") {
      p_moderate_severe <- input_params$p_moderate_severe_trt
    }  else {
      warning("Invalid treatment type (moderate_severe)") 
    }
    if (Trt == "Standard of care") {
      p_mild_severe <- input_params$p_mild_severe_soc
    } else if (Trt == "FX06 treatment") {
      p_mild_severe <- input_params$p_mild_severe_trt
    }  else {
      warning("Invalid treatment type (mild_severe)") 
    }
    if (Trt == "Standard of care") {
      p_MV_ICU_mild <- input_params$p_MV_ICU_mild_soc
    } else if (Trt == "FX06 treatment") {
      p_MV_ICU_mild <- input_params$p_MV_ICU_mild_trt
    } else {
      warning("Invalid treatment type (MV_ICU_mild)") 
    }
    if (Trt == "Standard of care") {
      p_MV_ICU_moderate <- input_params$p_MV_ICU_moderate_soc
    } else if (Trt == "FX06 treatment") {
      p_MV_ICU_moderate <- input_params$p_MV_ICU_moderate_trt
    } else { 
      warning("Invalid treatment type (MV_ICU_moderate") 
    }
    if (Trt == "Standard of care") {
      p_MV_ICU_severe <- input_params$p_MV_ICU_severe_soc
    } else if (Trt == "FX06 treatment") {
      p_MV_ICU_severe <- input_params$p_MV_ICU_severe_trt
    } else {
      warning("Invalid treatment type (MV_ICU_severe)") 
    }
    
    # create matrix of state transition probabilities
    m_p_t           <- matrix(0, nrow = input_params$n_states, ncol = input_params$n_i)  
    # give the state names to the rows
    rownames(m_p_t) <- input_params$v_names_states                               
    
    hr_RTW_all       <- inner_join(df_X, df_hr_iculos_rtw, by = "n_cycles_ICU")
    hr_RTW           <- hr_RTW_all[M_t == "REH", "hr_RTW"]
    
    p_REH_D_all      <- inner_join(df_X, df_p_REH_D, by = "age_group_id")
    p_REH_D          <- p_REH_D_all[M_t == "REH", "p_REH_D"]
    
    p_REC_D_all      <- inner_join(df_X, df_p_REC_D, by = "Age")
    p_REC_D          <- p_REC_D_all[M_t == "REC", "p_mort_d"]
    
    #### mild #### 
    # transition probabilities when hospitalised in ICU on MV
    m_p_t["MV_mild",      M_t == "MV_mild"]  <- (1 - input_params$p_D_mild - max(p_MV_ICU_mild[df_X$n_cycles_MV]) - max(p_mild_moderate[df_X$mild_moderate_YN]) - max(p_mild_severe[df_X$mild_severe_YN]))
    m_p_t["ICU_mild",     M_t == "MV_mild"]  <- max(p_MV_ICU_mild[df_X$n_cycles_MV])
    m_p_t["D",            M_t == "MV_mild"]  <- input_params$p_D_mild    
    m_p_t["MV_moderate",  M_t == "MV_mild"] <- max(p_mild_moderate[df_X$mild_moderate_YN])
    m_p_t["MV_severe",    M_t == "MV_mild"] <- max(p_mild_severe[df_X$mild_severe_YN])
    
    # transition probabilities when hospitalised in ICU without MV
    m_p_t["ICU_mild",     M_t == "ICU_mild"] <- 1 - input_params$p_D_mild - input_params$p_ICU_GW_mild - max(p_mild_moderate[df_X$mild_moderate_YN]) - max(p_mild_severe[df_X$mild_severe_YN])
    m_p_t["GW_mild",      M_t == "ICU_mild"] <- input_params$p_ICU_GW_mild        
    m_p_t["D",            M_t == "ICU_mild"] <- input_params$p_D_mild    
    m_p_t["ICU_moderate", M_t == "ICU_mild"] <- max(p_mild_moderate[df_X$mild_moderate_YN])
    m_p_t["ICU_severe",   M_t == "ICU_mild"] <- max(p_mild_severe[df_X$mild_severe_YN])
    
    # transition probabilities when hospitalised outside ICU
    m_p_t["GW_mild",      M_t == "GW_mild"] <- (1 - input_params$p_D_mild - input_params$p_GW_REH_mild)
    m_p_t["REH",          M_t == "GW_mild"] <- input_params$p_GW_REH_mild
    m_p_t["D",            M_t == "GW_mild"] <- input_params$p_D_mild    
    
    #### moderate ####
    # transition probabilities when hospitalised in ICU on MV
    m_p_t["MV_moderate",  M_t == "MV_moderate"]  <- (1 - input_params$p_D_moderate - max(p_MV_ICU_moderate[df_X$n_cycles_MV]) - max(p_moderate_severe[df_X$moderate_severe_YN]))
    m_p_t["ICU_moderate", M_t == "MV_moderate"]  <- max(p_MV_ICU_moderate[df_X$n_cycles_MV])
    m_p_t["D",            M_t == "MV_moderate"]  <- input_params$p_D_moderate 
    m_p_t["MV_severe",    M_t == "MV_moderate"] <- max(p_moderate_severe[df_X$moderate_severe_YN])
    
    # transition probabilities when hospitalised in ICU without MV
    m_p_t["ICU_moderate", M_t == "ICU_moderate"] <- (1 - input_params$p_D_moderate - input_params$p_ICU_GW_moderate - max(p_moderate_severe[df_X$moderate_severe_YN]))
    m_p_t["GW_moderate",  M_t == "ICU_moderate"] <- input_params$p_ICU_GW_moderate 
    m_p_t["D",            M_t == "ICU_moderate"] <- input_params$p_D_moderate    
    m_p_t["ICU_severe",   M_t == "ICU_moderate"] <- max(p_moderate_severe[df_X$moderate_severe_YN])
    
    # transition probabilities when hospitalised outside ICU
    m_p_t["GW_moderate",  M_t == "GW_moderate"] <- (1 - input_params$p_D_moderate - input_params$p_GW_REH_moderate)
    m_p_t["REH",          M_t == "GW_moderate"] <- input_params$p_GW_REH_moderate
    m_p_t["D",            M_t == "GW_moderate"] <- input_params$p_D_moderate
    
    #### severe ####
    # transition probabilities when hospitalised in ICU on MV
    m_p_t["MV_severe",    M_t == "MV_severe"]  <- (1 - input_params$p_D_severe - max(p_MV_ICU_severe[df_X$n_cycles_MV]))
    m_p_t["ICU_severe",   M_t == "MV_severe"]  <- max(p_MV_ICU_severe[df_X$n_cycles_MV])
    m_p_t["D",            M_t == "MV_severe"]  <- input_params$p_D_severe    
    
    # transition probabilities when hospitalised in ICU without MV
    m_p_t["ICU_severe",   M_t == "ICU_severe"] <- (1 - input_params$p_D_severe - input_params$p_ICU_GW_severe) 
    m_p_t["GW_severe",    M_t == "ICU_severe"] <- input_params$p_ICU_GW_severe 
    m_p_t["D",            M_t == "ICU_severe"] <- input_params$p_D_severe    
    
    # transition probabilities when hospitalised outside ICU
    m_p_t["GW_severe",    M_t == "GW_severe"] <- (1 - input_params$p_D_severe - input_params$p_GW_REH_severe) 
    m_p_t["REH",          M_t == "GW_severe"] <- input_params$p_GW_REH_severe
    m_p_t["D",            M_t == "GW_severe"] <- input_params$p_D_severe    
    
    
    # transition probabilities when in rehabilitation
    m_p_t["REH",          M_t == "REH"] <- 1 - p_REH_D - (input_params$p_REH_REC * hr_RTW)
    m_p_t["REC",          M_t == "REH"] <- input_params$p_REH_REC * hr_RTW
    m_p_t["D",            M_t == "REH"] <- p_REH_D
    
    # transition probabilities when recovered
    m_p_t["REC",          M_t == "REC"] <- 1 - p_REC_D
    m_p_t["D",            M_t == "REC"] <- p_REC_D 
    
    # transition probabilities when dead
    m_p_t["D",            M_t == "D"]  <- 1 
    
    return(t(m_p_t))
  }
  
  Costs <- function (M_t, Trt, input_params, df_X, df_c) {
    # Arguments:
    # M_t: health state occupied at cycle t (character variable)
    # Returns: 
    # costs accrued in this cycle
    # Trt:  treatment
    
    # Treatment specific transition costs
    if (Trt == "Standard of care") {
      c_trt <- 0
    } else if (Trt == "FX06 treatment") {
      c_trt <- df_c$value[which(df_c$type == "c_trt")]
    } 
    
    c_t <- c()
    
    c_t[M_t %in% c("MV_mild", "MV_moderate", "MV_severe")] <- if_else(Trt == "FX06 treatment" & df_X$trt_YN == 1, df_c$value[which(df_c$type == "c_MV")] + df_c$value[which(df_c$type == "c_trt")], df_c$value[which(df_c$type == "c_MV")])
    c_t[M_t %in% c("ICU_mild", "ICU_moderate", "ICU_severe")] <- if_else(Trt == "FX06 treatment" & df_X$trt_YN == 1, df_c$value[which(df_c$type == "c_ICU")] + df_c$value[which(df_c$type == "c_trt")], df_c$value[which(df_c$type == "c_ICU")])
    c_t[M_t %in% c("GW_mild", "GW_moderate", "GW_severe")] <- if_else(Trt == "FX06 treatment" & df_X$trt_YN == 1, df_c$value[which(df_c$type == "c_GW")] + df_c$value[which(df_c$type == "c_trt")], df_c$value[which(df_c$type == "c_GW")])
    c_t[M_t == "REH"]  <- df_c$value[which(df_c$type == "c_REH")]
    c_t[M_t == "REC"]  <- df_c$value[which(df_c$type == "c_REC")]
    c_t[M_t == "D"]    <- df_c$value[which(df_c$type == "c_D")]
    
    return(c_t)  # return costs accrued this cycle
  }
  
  MicroSim <- function(Trt, seed, input_params, df_X, df_c, df_hr_iculos_rtw, df_p_REH_D, df_p_REC_D){
    
    set.seed(seed)
    n_states <- length(input_params$v_names_states) # the number of health states
    
    m_M <- m_C <-  matrix(nrow = input_params$n_i, ncol = input_params$n_cycles + 1, 
                          dimnames = list(paste("ind"  , 1:input_params$n_i, sep = " "), 
                                          paste("cycle", 0:input_params$n_cycles, sep = " ")))  
    
    m_M [, 1] <- as.character(input_params$df_X$M_init) # initial health state at cycle 0 for individual i
    # calculate costs per individual during cycle 0
    m_C[, 1]  <- Costs(m_M[, 1], Trt=Trt, input_params = input_params, df_X = df_X, df_c = df_c)     
    
    # open a loop for time running cycles 1 to n_cycles 
    for (t in 1:input_params$n_cycles) {
      # calculate the transition probabilities for the cycle based on health state t
      m_P <- Probs(m_M[, t], t, Trt, input_params, df_X, df_hr_iculos_rtw, df_p_REH_D, df_p_REC_D)
      # check if transition probabilities are between 0 and 1
      #check_transition_probability(m_P, verbose = F)
      # check if each of the rows of the transition probabilities matrix sum to one
      #check_sum_of_transition_array(m_P, n_rows = input_params$n_i, n_cycles = input_params$n_cycles, verbose = F)
      
      # sample the next health state and store that state in matrix m_M
      m_M[, t + 1]  <- samplev(m_P, 1)    
      # calculate costs per individual during cycle t + 1
      m_C[, t + 1]  <- Costs(m_M[, t + 1], Trt=Trt, input_params = input_params, df_X = df_X, df_c = df_c)     
      
      # update time since illness onset for t + 1 
      df_X$mild_moderate_YN <- df_X$mild_severe_YN <- if_else(t+1 < 7 & df_X$Severity == "mild", 1, 2)
      df_X$moderate_severe_YN <- if_else(t+1 < 7 & df_X$Severity == "moderate", 1, 2)
      
      df_X$n_cycles_MV <- if_else(grepl("MV", m_M[, t + 1]), df_X$n_cycles_MV +1, df_X$n_cycles_MV)
      
      df_X$n_cycles_ICU <- if_else(m_M[, t + 1] %in% c("MV_mild", "MV_moderate", "MV_severe", "ICU_mild", "ICU_moderate", "ICU_severe"), 
                                   df_X$n_cycles_ICU +1,
                                   df_X$n_cycles_ICU)
      
      df_X$trt_YN <- if_else(m_M[, t + 1] %in% c("MV_mild", "MV_moderate", "MV_severe",
                                                 "ICU_mild", "ICU_moderate", "ICU_severe",
                                                 "GW_mild", "GW_moderate", "GW_severe")
                             & t+1 < 5
                             , 1
                             , 0)
      
      # Display simulation progress
      if(t/(input_params$n_cycles/10) == round(t/(input_params$n_cycles/10), 0)) { # display progress every 10%
        cat('\r', paste(t/input_params$n_cycles * 100, "% done", sep = " "))
      }
      
    } # close the loop for the time points 
    
    tc      <- rowSums(m_C)   # total cost per individual
    tc_hat  <- mean(tc)       # average cost
    tc_tot  <- sum(tc)
    # store the results from the simulation in a list
    results <- list(df_X = df_X, m_M = m_M, m_C = m_C, tc = tc, tc_hat = tc_hat, tc_tot = tc_tot)   
    
    return(results)  # return the results
  }
  
  health_state_trace_plot <- function(Trt, m_trace, v_names_states, n_i){
    m_TR <- t(apply(m_trace, 2, function(x) table(factor(x, levels = v_names_states, ordered = TRUE))))
    # m_TR <- m_TR / n_i                                 # calculate the proportion of individuals
    #m_TR <- m_TR / nrow(outcomes_SoC$m_M)
    colnames(m_TR) <- v_names_states                  # name the rows of the matrix
    rownames(m_TR) <- 1:(ncol(m_trace)) # name the columns of the matrix
    # Plot trace of first health state
    
    df_m_TR <- m_TR %>% as.data.frame(m_TR) %>% mutate(day = row_number())
    pivot_m_TR <- df_m_TR %>% 
      pivot_longer(cols = c("MV_mild":"D"),
                   names_to = "Severity",
                   values_to = "Prop")
    
    pivot_m_healthstate <- pivot_m_TR %>% 
      mutate(healthstate = case_when(Severity %in% c("MV_mild", "MV_moderate", "MV_severe") ~ "MV",
                                     Severity %in% c("ICU_mild", "ICU_moderate", "ICU_severe") ~ "ICU",
                                     Severity %in% c("GW_mild", "GW_moderate", "GW_severe") ~ "GW",
                                     Severity == "REH" ~ "REH",
                                     Severity == "REC" ~ "REC",
                                     Severity == "D" ~ "D")
      )
    
    pivot_m_healthstate <- pivot_m_healthstate %>% 
      select(-Severity) %>% 
      group_by(day, healthstate) %>% 
      summarize(Proportion = sum(Prop)/n_i)
    
    plot <- pivot_m_healthstate %>% ggplot(aes(x=day,
                                               y=Proportion,
                                               group = healthstate,
                                               color = healthstate
    ), ylim = c(0, 1)) +
      geom_line(linewidth = 1) +
      labs(title = str_c("Health state trace for ", Trt)) 
    #+ transition_reveal(day)
  }
  
  costs_state_df <- function(m_M, m_C){
    
    df_c_TR <- as.data.frame(m_C) %>% mutate(Individual = row_number())
    df_m_TR <- as.data.frame(m_M) %>% mutate(Individual = row_number())
    
    pivot_c_TR <- df_c_TR %>% 
      pivot_longer(cols = c("cycle 0":"cycle 60"),
                   names_to = "Day",
                   values_to = "Costs")
    pivot_m_TR <- df_m_TR %>% 
      pivot_longer(cols = c("cycle 0":"cycle 60"),
                   names_to = "Day",
                   values_to = "State") %>% 
      mutate(State = case_when(State %in% c("MV_mild", "MV_moderate", "MV_severe") ~ "MV",
                               State %in% c("ICU_mild", "ICU_moderate", "ICU_severe") ~ "ICU",
                               State %in% c("GW_mild", "GW_moderate", "GW_severe") ~ "GW",
                               State == "REH" ~ "REH",
                               State == "REC" ~ "REC",
                               State == "D" ~ "D"))
    
    df_c_s  <- inner_join(pivot_c_TR, pivot_m_TR, by = c("Individual", "Day"))
    df_c_s  <- df_c_s %>% mutate(Day = as.numeric(str_extract(Day, "\\d+")))
  }
  
  avg_costs_df <- function(df_c_s){
    mean_c_s <- df_c_s %>%
      select(!Individual, !Day) %>%
      group_by(State) %>% 
      summarize(avg_costs = round(mean(Costs),2)) %>% 
      filter(avg_costs > 0)
  }
    
  tot_costs_df <- function(Trt, df_c_s){
    total_c_s <- df_c_s %>%
      select(!Individual) %>%
      group_by(State, Day) %>% 
      summarize(Costs_EUR = round(sum(Costs),2)) %>% 
      filter(Costs_EUR > 0)
    
    #plot <- total_c_s %>% ggplot(aes(x=Day,
    #           y=Costs_EUR,
    #           group = State,
    #           color = State)) +
    #  geom_line(linewidth = 1) +
    #  xlab("Days") +
    #  ylab("Costs in Euros") +
    #  labs(title = str_c("Total costs per health state for ", Trt),
    #       colour = "Health state")
  } 
  
  run_simulation <- function(country_analysis, n_i, p_i_mild, p_i_moderate, r_i_MV, c_trt, c_MV, c_ICU, c_GW, c_REH, c_REC, c_D, hr_mild_moderate, hr_moderate_severe, hr_mild_severe, hr_MV_ICU_mild, hr_MV_ICU_moderate, hr_MV_ICU_severe){
    input_params <- input_parameters(country_analysis, n_i, p_i_mild, p_i_moderate, r_i_MV, c_trt, c_MV, c_ICU, c_GW, c_REH, c_REC, c_D, hr_mild_moderate, hr_moderate_severe, hr_mild_severe, hr_MV_ICU_mild, hr_MV_ICU_moderate, hr_MV_ICU_severe)
    df_hr_iculos_rtw <- as.data.frame(input_params$df_hr_iculos_rtw)
    df_p_REH_D <- as.data.frame(input_params$df_p_REH_D)
    df_p_REC_D <- as.data.frame(input_params$df_p_REC_D)
    df_X <- as.data.frame(input_params$df_X)
    df_c <- as.data.frame(input_params$df_c)
    outcomes_SoC   <- MicroSim(Trt="Standard of care", seed = 77, input_params = input_params, df_X = df_X, df_c = df_c, df_hr_iculos_rtw = df_hr_iculos_rtw, df_p_REH_D = df_p_REH_D, df_p_REC_D = df_p_REC_D) 
    outcomes_trt   <- MicroSim(Trt="FX06 treatment", seed = 77, input_params = input_params, df_X = df_X, df_c = df_c, df_hr_iculos_rtw = df_hr_iculos_rtw, df_p_REH_D = df_p_REH_D, df_p_REC_D = df_p_REC_D)
    df_c_s_soc     <- costs_state_df(outcomes_SoC$m_M, outcomes_SoC$m_C) 
    df_c_s_trt     <- costs_state_df(outcomes_trt$m_M, outcomes_trt$m_C) 
    traceplot_soc  <- health_state_trace_plot(Trt = "Standard of care", outcomes_SoC$m_M, v_names_states = input_params$v_names_states, n_i = input_params$n_i)
    traceplot_trt  <- health_state_trace_plot(Trt = "FX06 treatment", outcomes_trt$m_M, v_names_states = input_params$v_names_states, n_i = input_params$n_i)
    avg_costs_soc  <- avg_costs_df(df_c_s_soc) 
    avg_costs_trt  <- avg_costs_df(df_c_s_trt)
    tot_costs_soc  <- tot_costs_df(Trt = "Standard of care", df_c_s_soc) 
    tot_costs_trt  <- tot_costs_df(Trt = "FX06 treatment", df_c_s_trt)
    
    outcomes <- list(outcomes_SoC = outcomes_SoC,
                     outcomes_trt = outcomes_trt,
                     plot_trace_soc = traceplot_soc,
                     plot_trace_trt = traceplot_trt,
                     avg_costs_soc = avg_costs_soc,
                     avg_costs_trt = avg_costs_trt,
                     tot_costs_soc = tot_costs_soc,
                     tot_costs_trt = tot_costs_trt)
    return(outcomes)
  }
  
  outcomes_st <- run_simulation(country_analysis = country_analysis, n_i = n_i, p_i_mild = p_i_mild, p_i_moderate = p_i_moderate, r_i_MV = r_i_MV, c_trt = c_trt, c_MV = c_MV, c_ICU = c_ICU, c_GW = c_GW, c_REH = c_REH, c_REC = c_REC, c_D = c_D, hr_mild_moderate = hr_mild_moderate, hr_moderate_severe = hr_moderate_severe, hr_mild_severe = hr_mild_severe, hr_MV_ICU_mild = hr_MV_ICU_mild, hr_MV_ICU_moderate = hr_MV_ICU_moderate, hr_MV_ICU_severe = hr_MV_ICU_severe)
  return(outcomes_st)
}
