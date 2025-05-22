# Define the path to the environment image
environment_image_path <- "environment_image.RData"
source("src/global.R")
# Check if environment image exists
if (file.exists(environment_image_path)) {
  # Load the environment image if it exists
  message("Loading pre-trained model from saved environment...")
  #load(environment_image_path)
  message("Environment loaded successfully!")
} else {
  # Run the full pipeline if environment image doesn't exist
  message("No saved environment found. Running full data processing pipeline...")
  message("[1/4]: Loading Datasets")
  source("src/load_data.R")
  message("[2/4]: Performing EDA")
  source("src/eda.R")
  message("[3/4]: Building Model")
  source("src/model.R")
  
  # Save the environment for future use
  message("Saving environment for faster startup in the future...")
  save.image(environment_image_path)
  message("Environment saved successfully!")
}

# Load UI and server components (always needed regardless of environment loading)
message("[4/4]: Loading UI and server")
source("src/server/server.R")
source("src/ui/ui.R")

# Launch the Shiny app
shinyApp(
  ui = ui, 
  server = server
)