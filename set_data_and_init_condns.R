set_data_and_init_condns <- function(str){
  if(str == "SF"){
    source("set_SF_data_and_init_condns.R")
  } else if (str=="Seattle"){
    source("set_Seattle_data_and_init_condns.R")
  } else if (str=="Seattle_A"){
    source("set_Seattle_A_data_and_init_condns.R")
  } else if (str=="Seattle_B"){
    source("set_Seattle_B_data_and_init_condns.R")
  } else if (str=="Seattle_C"){
    source("set_Seattle_C_data_and_init_condns.R")
  } else if (str=="Boston"){
    source("set_Boston_data_and_init_condns.R")
  }
} 