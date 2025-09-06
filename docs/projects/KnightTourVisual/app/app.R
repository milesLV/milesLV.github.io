library(tidyverse)
library(Matrix)
library(purrr)
library(expm)
library(grid)
library(ggplot2)
library(stringr)
library(shiny)
library(rsconnect)

'%notin%' <- Negate('%in%')


generate_combinations <- function(ops, nums) {
  # Convert to numeric if character
  nums <- as.numeric(nums)
  
  # Generate all permutations of numbers
  permute <- function(v) {
    if (length(v) == 1) return(list(v))
    out <- list()
    for (i in seq_along(v)) {
      rest <- v[-i]
      for (p in permute(rest)) {
        out <- append(out, list(c(v[i], p)))
      }
    }
    return(out)
  }
  
  # Parse operations into functions
  parse_op <- function(op_string) {
    if (op_string == "+") {
      return(function(x) x)
    } else if (op_string == "-") {
      return(function(x) -x)
    } else if (grepl("^\\^", op_string)) {
      n <- as.numeric(sub("^\\^", "", op_string))
      return(function(x) x^n)
    } else if (grepl("^/", op_string)) {
      n <- as.numeric(sub("^/", "", op_string))
      return(function(x) x / n)
    } else {
      stop(paste("Unknown operation:", op_string))
    }
  }
  
  # Convert all ops to functions
  op_funcs <- lapply(ops, parse_op)
  
  # Get all combinations of operations (each element gets an op)
  op_combos <- expand.grid(rep(list(seq_along(op_funcs)), length(nums)))
  
  # Generate number permutations
  num_perms <- permute(nums)
  
  result <- list()
  
  for (perm in num_perms) {
    for (i in 1:nrow(op_combos)) {
      indices <- as.integer(op_combos[i, ])
      funcs <- op_funcs[indices]
      transformed <- mapply(function(f, x) f(x), funcs, perm)
      result <- append(result, list(transformed))
    }
  }
  
  return(result)
}

# generate_combinations(c("+", "-"), c(1, 2))
makeKnightMoves <- function(position, board_dims = c(8,8)) {
  position <- as.numeric(position)
  numRows <- board_dims[1]
  numCols <- board_dims[2]
  
  ifValid <- function(placement) {
    rowInBounds <- placement[1] >= 1 && placement[1] <= numRows
    colInBounds <- placement[2] >= 1 && placement[2] <= numCols
    return(rowInBounds && colInBounds)
  }
  
  # Validate input
  if (!ifValid(position)) {
    stop("Please enter values within board dimensions")
  }
  
  movesList <- lapply(generate_combinations(c("+", "-"), c(1, 2)), function(a){
    newPosition <- c(position[1] + a[1], position[2] + a[2])
    if (ifValid(newPosition)) {
      return(newPosition)
    } else {
      return(NULL)
    }
  })
  
  return(compact(movesList))
}

makeKnightMoves(c(8,7), c(8,7))

#1\space\space\space\space 2\space\space\space\space3\space\space\space\space4\\
makeLatticeAdjMatrix <- function(numRow = 8, numCol = 8) {
  cellTotal <- numRow * numCol
  A <- matrix(0, nrow = cellTotal, ncol = cellTotal)
  index_from_coords <- function(row, col) {
    return((row - 1) * numCol + col)
  }
  
  for (rowIndex in 1:numRow) {
    for (colIndex in 1:numCol) {
      from_index <- index_from_coords(rowIndex, colIndex)
      moves <- makeKnightMoves(c(rowIndex, colIndex), board_dims = c(numRow, numCol))
      for (move in moves) {
        moveRow <- move[1]
        moveCol <- move[2]
        
        if (moveRow <= numRow && moveCol <= numCol) { # if it's in the board
          to_index <- index_from_coords(moveRow, moveCol)
          A[from_index, to_index] <- 1
          A[to_index, from_index] <- 1  # undirected
        }
      }
    }
  }
  
  # return(as(A, "sparseMatrix"))
  return(A)
}

plotKnightFromMatrix <- function(position = c(1, 1), nrows = 8, ncols = 8, exponent = 1, showPathsToExponent = FALSE, probability = FALSE) {
  
  if (length(position) != 2 || 
      !all(position >= 1) || 
      position[1] > nrows || 
      position[2] > ncols) {
    stop("Position must be within the board dimensions")
  }
  
  A <- makeLatticeAdjMatrix(nrows, ncols)
  
  if(showPathsToExponent && exponent != 1) {
    ATotal <- A
    for (i in 2:exponent) {
      ATotal <- ATotal + A %^% i
      
    }
    A <- ATotal
  } else {
    A <- A %^% exponent # exponentiating
  }
  
  # Convert origin to matrix index
  origin_index <- (position[1] - 1) * ncols + position[2]
  path_counts <- A[origin_index, ]
  
  # Map index → coordinates
  grid <- expand.grid(x = 1:ncols, y = 1:nrows)
  grid$index <- (grid$y - 1) * ncols + grid$x
  grid$paths <- path_counts[grid$index]
  
  cutOff = 5
  if (probability) {
    grid <- grid %>%
      mutate(totalPaths = sum(paths),
             paths = paths / totalPaths)
  }
  grid <- grid %>%
    mutate(numDigits = ifelse(str_detect(paths, "^0\\."), 
                              nchar(paths) - 3,
                              nchar(paths)),
           pathsFormatted = ifelse(numDigits %notin% c(0:cutOff),
                                   format(signif(paths, 3), scientific = TRUE), 
                                   paths))
  
  # Label the origin
  grid$Type <- ifelse(grid$index == origin_index, "Origin", "Other")
  origin_tile <- grid %>% filter(Type == "Origin")
  
  knightImage <- png::readPNG("knightImage.png")
  knightGrob <- rasterGrob(knightImage, interpolate = TRUE)
  
  readable = 16
  # Plot
  ggplot(grid, aes(x = x, y = y, fill = paths)) +
    geom_tile(color = "grey60") +
    annotation_custom(
      knightGrob,
      xmin = origin_tile$x - 0.5,
      xmax = origin_tile$x + 0.5,
      ymin = origin_tile$y - 0.5,
      ymax = origin_tile$y + 0.5
    ) +
    scale_fill_gradient(name = "# of Paths", low = "white", high = "firebrick", na.value = "white") +
    scale_x_continuous(breaks = 1:ncols, minor_breaks = NULL) + 
    scale_y_continuous(breaks = 1:nrows, minor_breaks = NULL) + 
    geom_text(aes(label = ifelse(paths > 0, pathsFormatted, "")), 
              color = "black", 
              size = ifelse(grid$numDigits <= 3, 
                            5, 
                            ifelse(grid$numDigits > cutOff, 
                                   4, 
                                   4.5))
    ) +
    coord_fixed() +
    labs(title = str_c("Possible Knight Moves after ", exponent," ", ifelse(exponent == 1, "Jump", "Consecutive Jumps"))) +
    guides(fill = guide_colorbar(barwidth = 25, barheight = 3)) + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.ticks.length = unit(0, "pt"),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 16, margin = margin(t = -5)),
          axis.text.y = element_text(size = 16, margin = margin(r = -5)),
          plot.title = element_text(size = readable+5, margin = margin(b = -10)),
          legend.text = element_text(size = readable-1),
          legend.title = element_text(size = readable+2)
    )
}
# plotKnightFromMatrix(c(6,5), 8, 6, exponent = 1, showPathsToExponent = TRUE, probability = TRUE)

ui <- fluidPage(
  sidebarLayout( # all of the input bars
    sidebarPanel(
      div(
        style = "font-size: 24px; font-weight: bold; margin-top: 0px;",
        textOutput("textBox")
      ),
      sliderInput("exponent", # id
                  label = "Choose # of Consecutive Jumps",
                  min = 1,
                  max = 25,
                  value = 1,
                  step = 1),
      sliderInput("rows",
                  label = "# of Rows",
                  min = 2,
                  max = 8,
                  value = 8,
                  step = 1),
      sliderInput("cols",
                  label = "# of Columns",
                  min = 2,
                  max = 8,
                  value = 8,
                  step = 1),
      radioButtons("showPaths",
                   label = "Show Pathing Squares?",
                   choices = list("Yes" = TRUE, "No" = FALSE),
                   selected = TRUE),
      radioButtons("mode", # make probability of count
                   label = "Number Mode",
                   choices = list("Count" = FALSE, "Probability" = TRUE),
                   selected = FALSE)
    ),
    mainPanel(
      plotOutput("knightPlot", click = "plot_click", width = "100%")
    )
  )
)

# Server logic to receive inputs
server <- function(input, output, session) {
  pos <- reactiveVal(c(5, 5))
  text <- reactiveVal("Click any square to move the knight")
  
  # === Reactive value to track last activity time ===
  lastActivity <- reactiveVal(Sys.time())
  
  # === Track whether auto animation is active ===
  autoAnimate <- reactiveVal(FALSE)
  
  # === Timestamp for programmatic changes to input$exponent ===
  lastProgrammatic <- reactiveVal(NULL)  # store Sys.time() for programmatic updates
  progIgnoreWindow <- 0.6                # seconds: window to treat exponent changes as programmatic
  
  # === Observe non-exponent inputs for activity (these are always user interactions) ===
  observeEvent(
    list(input$plot_click, input$rows, input$cols, input$showPaths, input$mode),
    {
      lastActivity(Sys.time())   # Reset activity
      autoAnimate(FALSE)         # Stop auto animation immediately
    }
  )
  
  # === Observe exponent changes but ignore ones that are very close to our programmatic updates ===
  observeEvent(input$exponent, {
    lp <- lastProgrammatic()
    if (!is.null(lp) && difftime(Sys.time(), lp, units = "secs") < progIgnoreWindow) {
      # This change was likely caused by our code -> ignore for "user activity"
      return()
    }
    # Otherwise treat as real user activity
    lastActivity(Sys.time())
    autoAnimate(FALSE)
  })
  
  # === Check inactivity every 500ms and start autoAnimate when >= 5s of inactivity ===
  observe({
    invalidateLater(500, session)
    isolate({
      if (difftime(Sys.time(), lastActivity(), units = "secs") >= 5 && !autoAnimate()) {
        autoAnimate(TRUE)
        
        # Immediately perform first increment and mark it programmatic so it won't count as activity
        lastProgrammatic(Sys.time())
        current <- isolate(input$exponent)
        if (current < 25) {
          updateSliderInput(session, "exponent", value = current + 1)
        } else {
          updateSliderInput(session, "exponent", value = 1)
        }
      }
    })
  })
  
  # === Increment exponent every 1.5 seconds while animating ===
  observe({
    req(autoAnimate())
    invalidateLater(1500, session)  # 1.5 second between increments
    
    isolate({
      # Before programmatic update, set timestamp so the activity observer ignores it
      lastProgrammatic(Sys.time())
      current <- input$exponent
      if (current < 25) {
        updateSliderInput(session, "exponent", value = current + 1)
      } else {
        updateSliderInput(session, "exponent", value = 1)
      }
    })
  })
  
  # === Update knight position on click ===
  observeEvent(input$plot_click, {
    click <- input$plot_click
    x <- max(min(round(click$x), input$cols), 1)
    y <- max(min(round(click$y), input$rows), 1)
    
    pos(c(y, x))
    text("")
  })
  
  # === Keep position within bounds ===
  observeEvent(input$rows, {
    current_pos <- pos()
    if (current_pos[1] > input$rows) pos(c(input$rows, current_pos[2]))
  })
  
  observeEvent(input$cols, {
    current_pos <- pos()
    if (current_pos[2] > input$cols) pos(c(current_pos[1], input$cols))
  })
  
  # === Render plot and text ===
  output$textBox <- renderText({ text() })
  
  output$knightPlot <- renderPlot({
    plotKnightFromMatrix(position = pos(), 
                         nrows = input$rows, 
                         ncols = input$cols, 
                         exponent = input$exponent, 
                         showPathsToExponent = as.logical(input$showPaths),
                         probability = as.logical(input$mode))
  }, height = 650, width = 650)
}

# Run the application
shinyApp(ui = ui, server = server)