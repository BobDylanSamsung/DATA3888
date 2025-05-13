source("src/global.R")
message("[1/4]: Loading Datasets")
source("src/load_data.R")
message("[2/4]: Perfoming EDA")
source("src/eda.R")
message("[3/4]: Building Model")
source("src/model.R")
message("[4/4]: Loading ui and server")
source("src/server.R")
source("src/ui.R")
shinyApp(
  ui = ui, 
  server = server, 
)
