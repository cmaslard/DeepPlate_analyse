# Function to evaluate color contrast
evaluate_contrast <- function(color) {
  # Convert color to LAB color space
  color_c_ratio <- colorspace::contrast_ratio(color,"white")
  
  # Return "black" if luminance is below a threshold, otherwise "white"
  if (color_c_ratio < 4.5) {
    return("black")
  } else {
    return("white")
  }
}



