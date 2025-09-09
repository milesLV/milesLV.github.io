library(dplyr)
library(Matrix)
library(grid)
library(ggplot2)
library(shiny)
library(shinylive)

'%notin%' <- Negate('%in%')

generate_combinations <- function(ops, nums) {
  nums <- as.numeric(nums)

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
  
  op_funcs <- lapply(ops, parse_op)
  op_combos <- expand.grid(rep(list(seq_along(op_funcs)), length(nums)))
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

makeKnightMoves <- function(position, board_dims = c(8,8)) {
  position <- as.numeric(position)
  numRows <- board_dims[1]
  numCols <- board_dims[2]
  
  ifValid <- function(placement) {
    rowInBounds <- placement[1] >= 1 && placement[1] <= numRows
    colInBounds <- placement[2] >= 1 && placement[2] <= numCols
    return(rowInBounds && colInBounds)
  }
  
  if (!ifValid(position)) {
    stop("Please enter values within board dimensions")
  }
  
  movesList <- lapply(generate_combinations(c("+", "-"), c(1, 2)), function(a){
    newPosition <- c(position[1] + a[1], position[2] + a[2])
    if (ifValid(newPosition)) {
      return(newPosition)
    } else {
      return(invisible())
    }
  })
  movesList[sapply(movesList, is.null)] <- NULL
  return(movesList)
}

makeKnightMoves(c(8,7), c(8,7))

matrix_power <- function(A, n) {
  if (!is.matrix(A) || ncol(A) != nrow(A)) {
    stop("A must be a square matrix")
  }
  if (!is.numeric(n) || length(n) != 1 || floor(n) != n) {
    stop("Exponent n must be a single integer")
  }
  
  if (n == 0) return(diag(nrow(A)))
  if (n < 0) {
    A <- solve(A)
    n <- -n
  }
  result <- diag(nrow(A))
  base <- A
  while (n > 0) {
    if (n %% 2 == 1) {
      result <- result %*% base
    }
    base <- base %*% base
    n <- n %/% 2
  }
  result
}

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
        
        if(moveRow <= numRow && moveCol <= numCol) {
          to_index <- index_from_coords(moveRow, moveCol)
          A[from_index, to_index] <- 1
          A[to_index, from_index] <- 1
        }
      }
    }
  }
  
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
      ATotal <- ATotal + matrix_power(A, i)
      
    }
    A <- ATotal
  } else {
    A <- matrix_power(A, exponent)
  }
  
  origin_index <- (position[1] - 1) * ncols + position[2]
  path_counts <- A[origin_index, ]
  
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
    mutate(numDigits = ifelse(grepl("^0\\.", paths),
                              nchar(paths) - 3,
                              nchar(paths)),
           pathsFormatted = ifelse(numDigits %notin% c(0:cutOff),
                                   format(signif(paths, 3), scientific = TRUE), 
                                   paths))
  
  grid$Type <- ifelse(grid$index == origin_index, "Origin", "Other")
  origin_tile <- grid %>% filter(Type == "Origin")
  
  knightImage <- png::readPNG("knightImage.png")
  knightGrob <- rasterGrob(knightImage, interpolate = TRUE)
  
  readable = 16
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
    labs(title = paste0("Possible Knight Moves after ", exponent," ", ifelse(exponent == 1, "Jump", "Consecutive Jumps"))) +
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

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      div(
        style = "font-size: 24px; font-weight: bold; margin-top: 0px;",
        textOutput("textBox")
      ),
      sliderInput("exponent",
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
      radioButtons("mode",
                   label = "Number Mode",
                   choices = list("Count" = FALSE, "Probability" = TRUE),
                   selected = FALSE),
      actionButton("toggleAnimate", "Toggle Auto Animate", class = "btn-primary")
    ),
    mainPanel(
      plotOutput("knightPlot", click = "plot_click", width = "100%")
    )
  )
)

server <- function(input, output, session) {
  pos <- reactiveVal(c(5, 5))
  text <- reactiveVal("Click any square to move the knight")
  autoAnimate <- reactiveVal(FALSE)
  lastProgrammatic <- reactiveVal(NULL)
  
  observeEvent(input$toggleAnimate, {
    autoAnimate(!autoAnimate())
  })
  
  observe({
    req(autoAnimate())
    invalidateLater(1500, session)
    
    isolate({
      lastProgrammatic(Sys.time())
      current <- input$exponent
      if (current < 25) {
        updateSliderInput(session, "exponent", value = current + 1)
      } else {
        updateSliderInput(session, "exponent", value = 1)
      }
    })
  })
  
  observeEvent(input$plot_click, {
    click <- input$plot_click
    x <- max(min(round(click$x), input$cols), 1)
    y <- max(min(round(click$y), input$rows), 1)
    
    pos(c(y, x))
    text("")
  })
  
  observeEvent(input$rows, {
    current_pos <- pos()
    if (current_pos[1] > input$rows) pos(c(input$rows, current_pos[2]))
  })
  
  observeEvent(input$cols, {
    current_pos <- pos()
    if (current_pos[2] > input$cols) pos(c(current_pos[1], input$cols))
  })
  
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

shinyApp(ui, server)